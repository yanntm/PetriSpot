/*
 * CtlDriver.h
 *
 * CTL properties from the command line: normalise each formula, decide the
 * state formulas at the initial marking, run the checker in rounds of
 * growing budgets under a wall clock per property, print the FORMULA lines
 * and, with --trace, the evidence trees.
 */
#ifndef PETRI_CLI_CTLDRIVER_H_
#define PETRI_CLI_CTLDRIVER_H_

#include <chrono>
#include <iostream>
#include <string>
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

/** Check the CTL properties; returns how many were decided. */
template<typename T>
  size_t runCtl (const Options &o, const SparsePetriNet<T> &pn, const std::vector<petri::expr::Property> &props)
  {
    using petri::ctl::Verdict;
    using clock = std::chrono::steady_clock;
    auto start = clock::now ();
    petri::walk::WalkNet<T> wnet (pn);
    const std::vector<std::string> &pnames = pn.getPnames ();
    const std::vector<std::string> &tnames = pn.getTnames ();
    size_t decided = 0;
    size_t open = props.size ();
    for (const auto &prop : props) {
      petri::expr::CtlFormula f = petri::expr::ctlNormalize (prop.ctl);
      if (!o.quiet) {
        std::cout << "CTL " << prop.name << " : ";
        f.print (std::cout, &pnames);
        std::cout << "  (" << f.size () << " nodes, depth " << f.depth () << ")" << std::endl;
      }
      petri::ctl::Cursor<T> init (wnet);
      if (f.isState ()) {
        bool holds = f.evalState (init.marking, init.deadlock ());
        std::cout << "FORMULA " << prop.name << " " << (holds ? "TRUE" : "FALSE")
            << " TECHNIQUES TOPOLOGICAL INITIAL_STATE" << std::endl;
        ++decided;
        --open;
        continue;
      }
      // the clock: an equal share of what is left of the total, or the timeout
      long long remainingMs = o.totalTime > 0
          ? o.totalTime * 1000LL - std::chrono::duration_cast<std::chrono::milliseconds> (clock::now () - start).count ()
          : static_cast<long long> (o.timeout) * 1000LL;
      long long shareMs = std::max<long long> (100, o.totalTime > 0 ? remainingMs / static_cast<long long> (open) : remainingMs);
      auto propStart = clock::now ();
      petri::ctl::Checker<T> checker (wnet, o.seed + 7919 * decided);
      petri::ctl::Budget budget;
      budget.runLength = o.ctlRunLength;
      budget.epsilon = o.epsilon;
      budget.saturate = o.ctlSaturate;
      budget.sample = o.sample ? o.sample : 8;
      budget.deadline = propStart + std::chrono::milliseconds (shareMs);
      Verdict v = Verdict::Unknown;
      petri::ctl::Evidence ev;
      unsigned round = 0;
      for (; round < o.ctlRounds && v == Verdict::Unknown && !checker.timedOut (); ++round) {
        // the work of a round grows tenfold: odd rounds widen the root search, even rounds the probe of a state
        uint64_t rootScale = 1, nestedScale = 1;
        for (unsigned r = 0; r < (round + 1) / 2; ++r) rootScale *= 10;
        for (unsigned r = 0; r < round / 2; ++r) nestedScale *= 10;
        budget.level = round;
        budget.huntSteps = o.ctlSteps * rootScale;
        budget.regionStates = o.ctlRegion * rootScale;
        // the probe of a state: a short hunt (a reachable operand is found in a few steps, an unreachable one
        // wastes the whole budget) and the region, which is where the proof is
        budget.nestedSteps = std::max<uint64_t> (100, o.ctlSteps / 10) * nestedScale;
        budget.nestedStates = o.ctlRegion * nestedScale;
        checker.setBudget (budget);
        ev = petri::ctl::Evidence ();
        v = checker.solve (init, f, &ev);
      }
      long ms = static_cast<long> (std::chrono::duration_cast<std::chrono::milliseconds> (clock::now () - propStart).count ());
      if (v != Verdict::Unknown) {
        std::cout << "FORMULA " << prop.name << " " << to_string (v) << " TECHNIQUES EXPLICIT CTL_WALK" << std::endl;
        ++decided;
        if (o.trace) {
          std::cout << "WITNESS " << prop.name << " (" << to_string (v) << ")\n";
          ev.print (std::cout, &tnames, 2);
          std::cout.flush ();
        }
      } else if (o.printUnknown) {
        std::cout << "UNKNOWN " << prop.name << std::endl;
      }
      if (!o.quiet) {
        std::cout << "CTL " << prop.name << " " << to_string (v) << " after " << round << " round(s), " << ms << " ms: ";
        checker.stats ().print (std::cout);
        std::cout << std::endl;
      }
      --open;
    }
    return decided;
  }

} // namespace petri::cli

#endif /* PETRI_CLI_CTLDRIVER_H_ */
