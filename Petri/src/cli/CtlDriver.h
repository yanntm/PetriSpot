/*
 * CtlDriver.h
 *
 * CTL properties from the command line: normalise each formula, decide the
 * state formulas at the initial marking, run the checker in rounds of
 * growing budgets under a wall clock per property, --threads properties at a
 * time, print the FORMULA lines and, with --trace, the evidence trees.
 */
#ifndef PETRI_CLI_CTLDRIVER_H_
#define PETRI_CLI_CTLDRIVER_H_

#include <algorithm>
#include <atomic>
#include <chrono>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "cli/Options.h"
#include "core/SparsePetriNet.h"
#include "ctl/Checker.h"
#include "ctl/Cursor.h"
#include "ctl/Evidence.h"
#include "expr/CtlSimplify.h"
#include "expr/Property.h"
#include "walk/WalkNet.h"

namespace petri::cli
{

/** The check of one CTL property: its normal form, the rounds, the lines to print. */
template<typename T>
  struct CtlOutcome
  {
    petri::ctl::Verdict verdict = petri::ctl::Verdict::Unknown;
    bool atInitial = false; // a state formula decided at the initial marking
    unsigned rounds = 0;
    long ms = 0;
    petri::ctl::Evidence evidence;
    petri::ctl::Stats stats;
  };

/** Run the rounds on one property under its clock; single-threaded. */
template<typename T>
  CtlOutcome<T> checkCtl (const Options &o, const petri::walk::WalkNet<T> &wnet, const petri::expr::CtlFormula &f,
                          uint64_t seed, long long clockMs)
  {
    using petri::ctl::Verdict;
    using clock = std::chrono::steady_clock;
    CtlOutcome<T> out;
    petri::ctl::Cursor<T> init (wnet);
    if (f.isState ()) {
      out.verdict = f.evalState (init.marking, init.deadlock ()) ? Verdict::True : Verdict::False;
      out.atInitial = true;
      return out;
    }
    auto propStart = clock::now ();
    petri::ctl::Checker<T> checker (wnet, seed);
    petri::ctl::Budget budget;
    budget.runLength = o.ctlRunLength;
    budget.epsilon = o.epsilon;
    budget.saturate = o.ctlSaturate;
    budget.sample = o.sample ? o.sample : 8;
    budget.deadline = propStart + std::chrono::milliseconds (clockMs);
    for (; out.rounds < o.ctlRounds && out.verdict == Verdict::Unknown && !checker.timedOut (); ++out.rounds) {
      // the work of a round grows tenfold: odd rounds widen the root search, even rounds the probe of a state
      uint64_t rootScale = 1, nestedScale = 1;
      for (unsigned r = 0; r < (out.rounds + 1) / 2; ++r) rootScale *= 10;
      for (unsigned r = 0; r < out.rounds / 2; ++r) nestedScale *= 10;
      budget.level = out.rounds;
      budget.huntSteps = o.ctlSteps * rootScale;
      budget.regionStates = o.ctlRegion * rootScale;
      // the probe of a state: a short hunt (a reachable operand is found in a few steps, an unreachable one
      // wastes the whole budget) and the region, which is where the proof is
      budget.nestedSteps = std::max<uint64_t> (100, o.ctlSteps / 10) * nestedScale;
      budget.nestedStates = o.ctlRegion * nestedScale;
      checker.setBudget (budget);
      out.evidence = petri::ctl::Evidence ();
      out.verdict = checker.solve (init, f, &out.evidence);
    }
    out.ms = static_cast<long> (std::chrono::duration_cast<std::chrono::milliseconds> (clock::now () - propStart).count ());
    out.stats = checker.stats ();
    return out;
  }

/**
 * Check the CTL properties, --threads of them at a time (each check is
 * single-threaded); returns how many were decided. The clock of a property is
 * -t, or the share of --totalTime of its wave of --threads properties.
 */
template<typename T>
  size_t runCtl (const Options &o, const SparsePetriNet<T> &pn, const std::vector<petri::expr::Property> &props)
  {
    using petri::ctl::Verdict;
    petri::walk::WalkNet<T> wnet (pn);
    const std::vector<std::string> &pnames = pn.getPnames ();
    const std::vector<std::string> &tnames = pn.getTnames ();
    std::vector<petri::expr::CtlFormula> forms;
    for (const auto &prop : props) forms.push_back (petri::expr::ctlNormalize (prop.ctl));
    if (!o.quiet) {
      for (size_t i = 0; i < props.size (); ++i) {
        std::cout << "CTL " << props[i].name << " : ";
        forms[i].print (std::cout, &pnames);
        std::cout << "  (" << forms[i].size () << " nodes, depth " << forms[i].depth () << ")" << std::endl;
      }
    }
    unsigned threads = std::max<unsigned> (1, o.threads);
    size_t waves = (props.size () + threads - 1) / threads;
    long long clockMs = o.totalTime > 0 ? std::max<long long> (100, o.totalTime * 1000LL / static_cast<long long> (waves))
                                        : static_cast<long long> (o.timeout) * 1000LL;
    std::mutex outMutex;
    std::atomic<size_t> next { 0 }, decided { 0 };
    auto worker = [&] () {
      for (size_t i = next++; i < props.size (); i = next++) {
        CtlOutcome<T> out = checkCtl (o, wnet, forms[i], o.seed + 7919 * static_cast<uint64_t> (i), clockMs);
        std::lock_guard<std::mutex> lock (outMutex);
        const std::string &name = props[i].name;
        if (out.verdict != Verdict::Unknown) {
          ++decided;
          std::cout << "FORMULA " << name << " " << to_string (out.verdict) << " TECHNIQUES "
              << (out.atInitial ? "TOPOLOGICAL INITIAL_STATE" : "EXPLICIT CTL_WALK") << std::endl;
          if (o.trace && !out.atInitial) {
            std::cout << "WITNESS " << name << " (" << to_string (out.verdict) << ")\n";
            out.evidence.print (std::cout, &tnames, 2);
            std::cout.flush ();
          }
        } else if (o.printUnknown) {
          std::cout << "UNKNOWN " << name << std::endl;
        }
        if (!o.quiet && !out.atInitial) {
          std::cout << "CTL " << name << " " << to_string (out.verdict) << " after " << out.rounds << " round(s), "
              << out.ms << " ms: ";
          out.stats.print (std::cout);
          std::cout << std::endl;
        }
      }
    };
    std::vector<std::thread> pool;
    for (unsigned k = 1; k < threads && k < props.size (); ++k) pool.emplace_back (worker);
    worker ();
    for (auto &th : pool) th.join ();
    return decided;
  }

} // namespace petri::cli

#endif /* PETRI_CLI_CTLDRIVER_H_ */
