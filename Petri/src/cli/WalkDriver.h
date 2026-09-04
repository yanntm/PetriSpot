/*
 * WalkDriver.h
 *
 * Reachability from the command line: load and print properties, build the
 * target set, run the sweep round and the focused rounds, print the FORMULA
 * lines as targets fall and the witnesses afterwards.
 */
#ifndef PETRI_CLI_WALKDRIVER_H_
#define PETRI_CLI_WALKDRIVER_H_

#include <algorithm>
#include <chrono>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "cli/Options.h"
#include "core/SparsePetriNet.h"
#include "expr/Property.h"
#include "expr/SexprPrinter.h"
#include "parse/PropertyFile.h"
#include "parse/sexpr/HintReader.h"
#include "walk/Portfolio.h"
#include "walk/SharedPool.h"
#include "walk/Target.h"
#include "walk/TargetSet.h"
#include "walk/WalkNet.h"

namespace petri::cli
{

/** Print a sparse marking with place names. */
template<typename T>
  void printMarking (std::ostream &os, const SparseArray<T> &m, const std::vector<std::string> &pnames)
  {
    for (size_t i = 0; i < m.size (); ++i) {
      if (i > 0) os << " ";
      os << pnames[m.keyAt (i)];
      if (m.valueAt (i) != 1) os << "=" << m.valueAt (i);
    }
  }

/** The strategy pool of the options: --strategies, else --strategy, else the default pool. */
inline std::vector<petri::walk::StrategySpec> strategyPool (const Options &o)
{
  std::string list = o.strategies;
  if (list.empty ()) {
    list = (o.strategy != "random" || o.threads == 1) ? o.strategy : "random,bestfirst,structural,relaxed";
  }
  return petri::walk::parseStrategySpecs (list, o.epsilon, o.stall, o.sample);
}

/** The MCC verdict line, printed the moment a target is claimed. */
template<typename T>
  void printFormula (const Options &o, const petri::walk::TargetSet<T> &targets, const petri::walk::Claim<T> &c)
  {
    std::cout << "FORMULA " << targets.name (c.target) << " " << targets.verdict (c.target) << " TECHNIQUES EXPLICIT "
        << (c.strategy == "random" ? "RANDOM_WALK" : "HEURISTIC_WALK") << (o.threads > 1 ? " PARALLEL_PROCESSING" : "")
        << std::endl;
  }

template<typename T>
  void printReports (const Options &o, const petri::walk::PortfolioResult<T> &res)
  {
    for (size_t i = 0; i < res.reports.size (); ++i) {
      const petri::walk::ThreadReport &rep = res.reports[i];
      const petri::walk::WalkStats &st = rep.stats;
      uint64_t ms = st.millis == 0 ? 1 : st.millis;
      if (o.quiet && !rep.found && res.found) continue;
      std::cout << "Thread " << i << " [" << rep.strategy << "] " << (rep.found ? "found a witness" : "stopped")
          << " after " << st.steps << " steps, " << st.resets << " resets (" << st.deadEnds << " dead ends, "
          << st.stalls << " stalls, " << st.poolRestarts << " from pool), " << st.millis << " ms ("
          << (st.steps / ms) << " steps/ms; " << (st.steps ? st.arcVisits / st.steps : 0) << " arc visits/step, "
          << (st.steps ? st.targetChecks / st.steps : 0) << " checks/step)";
      if (rep.strategy != "random") std::cout << ", best heuristic " << rep.minHeuristic;
      if (rep.claims > 1 || (rep.claims == 1 && !rep.found)) std::cout << ", " << rep.claims << " targets claimed";
      std::cout << "." << std::endl;
    }
  }

/** WITNESS lines (verified traces) or, when not quiet, the witness markings. */
template<typename T>
  void printWitnesses (const Options &o, const petri::walk::WalkNet<T> &wnet,
                       const petri::walk::TargetSet<T> &targets, const petri::walk::PortfolioResult<T> &res)
  {
    for (const auto &c : res.claims) {
      if (c.hasTrace) {
        const auto &tnames = wnet.getNet ().getTnames ();
        std::cout << "WITNESS " << targets.name (c.target) << " " << c.trace.size ();
        for (uint32_t t : c.trace) std::cout << " " << tnames[t];
        std::cout << std::endl;
      } else if (!o.quiet) {
        std::cout << "Witness marking for " << targets.name (c.target) << ": ";
        printMarking (std::cout, c.marking, wnet.getNet ().getPnames ());
        std::cout << std::endl;
      }
    }
  }

/** True iff every walker stopped because it exhausted the step budget. */
template<typename T>
  bool stepBound (const petri::walk::PortfolioResult<T> &res, const petri::walk::WalkBudget &budget)
  {
    if (budget.maxSteps == 0) return false;
    for (const auto &rep : res.reports)
      if (rep.stats.steps < budget.maxSteps) return false;
    return true;
  }

/**
 * Run the portfolio on targets toward focus (or as a sweep with NO_FOCUS);
 * FORMULA lines are printed as claims happen and the reports afterwards.
 */
template<typename T>
  petri::walk::PortfolioResult<T> runWalk (const Options &o, const petri::walk::WalkNet<T> &wnet, petri::walk::TargetSet<T> &targets,
                uint32_t focus, const petri::walk::WalkBudget &budget,
                const std::vector<petri::walk::StrategySpec> &specs)
  {
    std::unique_ptr<petri::walk::SharedPool<T>> pool;
    if (o.share > 0) pool = std::make_unique<petri::walk::SharedPool<T>> (o.share, o.shareProb);
    std::function<void (const petri::walk::Claim<T>&)> onClaim = [&] (const petri::walk::Claim<T> &c) {
      printFormula (o, targets, c);
    };
    petri::walk::PortfolioResult<T> res = petri::walk::runPortfolio (wnet, targets, focus, specs, o.threads, budget,
                                                                     o.seed, pool.get (), o.debugSteps, onClaim);
    printReports (o, res);
    if (pool) {
      std::cout << "Shared pool: " << pool->publishedCount () << " published, " << pool->drawnCount ()
          << " drawn, " << pool->evictedCount () << " evicted, " << pool->size () << " held." << std::endl;
    }
    printWitnesses (o, wnet, targets, res);
    if (focus != petri::walk::NO_FOCUS && !res.found) {
      std::cout << "No witness found for " << targets.name (focus) << "." << std::endl;
    }
    return res;
  }

/** Load the property file, keep --query if given. Throws std::string on error. */
template<typename T>
  std::vector<petri::expr::Property> loadProperties (const Options &o, const SparsePetriNet<T> &pn)
  {
    std::vector<petri::expr::Property> props = petri::loadPropertyFile (o.propsFile, pn,
                                                                          petri::propertySyntaxOf (o.propsSyntax));
    if (o.query >= 0) {
      if (static_cast<size_t> (o.query) >= props.size ()) {
        throw std::string ("--query=" + std::to_string (o.query) + " but the file holds "
            + std::to_string (props.size ()) + " properties.");
      }
      props = { props[static_cast<size_t> (o.query)] };
    }
    if (!o.hintsFile.empty ()) {
      size_t n = petri::sexpr::attachHints (props, petri::sexpr::loadHints (o.hintsFile, pn));
      std::cout << "Hints: " << n << " properties have a Parikh vector." << std::endl;
    }
    return props;
  }

/** --printProps: infix with the goal, or one s-expression form per line. */
template<typename T>
  void printProperties (const std::vector<petri::expr::Property> &props, const SparsePetriNet<T> &pn,
                        const std::string &format)
  {
    const std::vector<std::string> *pnames = format == "sexpr-index" ? nullptr : &pn.getPnames ();
    for (const auto &prop : props) {
      if (format == "infix") {
        prop.print (std::cout, pnames);
        std::cout << "\n  goal (" << prop.goal ().size () << " nodes) : ";
        prop.goal ().print (std::cout, pnames);
      } else {
        petri::expr::printSexpr (std::cout, prop, pnames);
      }
      std::cout << "\n";
    }
    std::cout.flush ();
  }

/** The running maxima of the bound targets (progress measure of a round). */
template<typename T>
  std::vector<long long> boundValues (const petri::walk::TargetSet<T> &targets)
  {
    std::vector<long long> v;
    for (uint32_t i = 0; i < targets.size (); ++i)
      if (targets.target (i).isBound ()) v.push_back (targets.bestValue (i));
    return v;
  }

/**
 * The target set of the supported, non-trivial properties. Unsupported ones
 * are reported, trivially false goals are answered here.
 */
template<typename T>
  petri::walk::TargetSet<T> makeTargets (const std::vector<petri::expr::Property> &props, const SparsePetriNet<T> &pn)
  {
    using petri::expr::PropertyKind;
    std::vector<petri::walk::Target<T>> targets;
    std::vector<std::string> names, verdicts;
    for (const auto &prop : props) {
      if (prop.kind == PropertyKind::Unsupported) {
        std::cout << "Skipping " << prop.name << " : " << prop.comment << std::endl;
        continue;
      }
      if (prop.kind == PropertyKind::Deadlock) {
        targets.push_back (petri::walk::Target<T>::deadlockTarget ());
      } else if (prop.kind == PropertyKind::Bound) {
        targets.push_back (petri::walk::Target<T>::boundTarget (prop.boundForm (), prop.boundHint));
      } else {
        petri::expr::Expression goal = prop.goal ();
        if (goal.kind == petri::expr::Expression::Kind::False) {
          std::cout << "FORMULA " << prop.name << " " << (prop.kind == PropertyKind::Invariant ? "TRUE" : "FALSE")
              << " TECHNIQUES TOPOLOGICAL TRIVIAL" << std::endl;
          continue;
        }
        targets.push_back (petri::walk::Target<T> (std::move (goal)));
      }
      targets.back ().setHint (prop.hint);
      names.push_back (prop.name);
      verdicts.push_back (prop.verdictIfReached ());
    }
    return petri::walk::TargetSet<T> (pn.getPlaceCount (), std::move (targets), std::move (names),
                                      std::move (verdicts));
  }

/**
 * Walk every property of the file: a sweep round (random walks, all targets
 * checked at once) when at least two are open, then focused rounds with a
 * per-property budget growing tenfold per round under --totalTime, or the -t
 * timeout for each otherwise. The rounds stop early when one of them solved
 * nothing, improved no bound, and every walk in it ended on the step budget:
 * more time would change nothing. Every bound target ends with a BOUND line.
 */
template<typename T>
  void runProperties (const Options &o, const SparsePetriNet<T> &pn)
  {
    using petri::walk::NO_FOCUS;
    std::vector<petri::expr::Property> props = loadProperties (o, pn);
    if (o.printProps) {
      printProperties (props, pn, o.printPropsFormat);
      return;
    }
    petri::walk::TargetSet<T> targets = makeTargets (props, pn);
    if (targets.size () == 0) return;
    petri::walk::WalkBudget budget = o.budget;
    budget.recordTrace = o.trace;
    petri::walk::WalkNet<T> wnet (pn);
    std::vector<petri::walk::StrategySpec> specs = strategyPool (o);
    std::vector<petri::walk::StrategySpec> hintSpecs = petri::walk::parseStrategySpecs (o.hintStrategies, o.epsilon,
                                                                                         o.stall, o.sample);
    std::vector<petri::walk::StrategySpec> randomSpecs = petri::walk::parseStrategySpecs ("random", 0, 0, 0);

    auto walkStart = std::chrono::steady_clock::now ();
    auto elapsedMs = [&] () { return millisSince (walkStart); };
    const long totalMs = o.totalTime * 1000;

    if (o.sweepTime > 0 && targets.openCount () >= 2) {
      long ms = o.sweepTime * 1000;
      if (totalMs > 0) ms = std::min (ms, totalMs);
      std::cout << "Sweep: " << targets.openCount () << " open properties, " << ms << " ms of random walks." << std::endl;
      budget.timeoutMillis = static_cast<uint64_t> (ms);
      runWalk (o, wnet, targets, NO_FOCUS, budget, randomSpecs);
    }

    long perProperty = (totalMs > 0 ? o.roundTime : o.timeout) * 1000; // milliseconds
    for (int round = 1; targets.openCount () > 0; ++round) {
      std::vector<uint32_t> open = targets.openTargets ();
      if (totalMs > 0) {
        long left = totalMs - elapsedMs ();
        // last round when the geometric budget would exceed what is left
        if (perProperty * static_cast<long> (open.size ()) >= left) perProperty = left / static_cast<long> (open.size ());
        if (perProperty < 20) break; // not worth a walk
        std::cout << "Round " << round << ": " << open.size () << " open properties, " << perProperty
            << " ms each, " << left << " ms left." << std::endl;
      }
      budget.timeoutMillis = static_cast<uint64_t> (perProperty);
      size_t openBefore = targets.openCount ();
      std::vector<long long> boundsBefore = boundValues (targets);
      bool allStepBound = true;
      for (uint32_t k : open) {
        if (targets.isSolved (k)) continue; // claimed on the way to another focus
        petri::walk::PortfolioResult<T> res = runWalk (o, wnet, targets, k, budget,
                                                       targets.target (k).hasHint () ? hintSpecs : specs);
        if (!stepBound (res, budget)) allStepBound = false;
        if (totalMs > 0 && elapsedMs () >= totalMs) break;
      }
      if (totalMs <= 0) break;
      if (allStepBound && targets.openCount () == openBefore && boundValues (targets) == boundsBefore) {
        std::cout << "Round " << round << " solved nothing and every walk ended on the step budget: stopping."
            << std::endl;
        break;
      }
      perProperty *= 10;
    }
    for (uint32_t i = 0; i < targets.size (); ++i) {
      if (targets.target (i).isBound ()) {
        std::cout << "BOUND " << targets.name (i) << " " << targets.bestValue (i) << std::endl;
      }
    }
    if (o.printUnknown) {
      for (uint32_t k : targets.openTargets ())
        if (!targets.target (k).isBound ()) std::cout << "UNKNOWN " << targets.name (k) << std::endl;
    }
  }

/** --findDeadlock: one deadlock target with the -t timeout. */
template<typename T>
  void runDeadlock (const Options &o, const SparsePetriNet<T> &pn)
  {
    petri::walk::WalkBudget budget = o.budget;
    budget.timeoutMillis = static_cast<uint64_t> (o.timeout) * 1000;
    budget.recordTrace = o.trace;
    petri::walk::WalkNet<T> wnet (pn);
    std::vector<petri::walk::StrategySpec> specs = petri::walk::parseStrategySpecs ("random", o.epsilon, o.stall, o.sample);
    std::vector<petri::walk::Target<T>> tg { petri::walk::Target<T>::deadlockTarget () };
    petri::walk::TargetSet<T> targets (pn.getPlaceCount (), std::move (tg), { "ReachabilityDeadlock" }, { "TRUE" });
    runWalk (o, wnet, targets, 0, budget, specs);
  }

} // namespace petri::cli

#endif /* PETRI_CLI_WALKDRIVER_H_ */
