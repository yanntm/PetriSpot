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
#include <limits>
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
#include "invariants/InvariantMiddle.h"
#include "walk/Components.h"
#include "walk/Knowledge.h"
#include "walk/Portfolio.h"
#include "walk/RestartPolicy.h"
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

/**
 * The processes of the net from its P-semiflows, when a strategy of the pool
 * wants them ("sync"), when the sweep asks for them, or when the sweep is
 * left to decide ("auto") with a few seconds of budget: the invariant engine
 * runs in-process, time-boxed, and the decomposition is printed. Null
 * otherwise.
 */
template<typename T>
  std::unique_ptr<petri::walk::Components<T>> buildComponents (const Options &o, const SparsePetriNet<T> &pn,
                                                               const petri::walk::WalkNet<T> &wnet,
                                                               const std::vector<petri::walk::StrategySpec> &specs,
                                                               size_t openTargets)
  {
    bool wanted = o.sweepChoice == "sync"
        || (o.sweepChoice == "auto" && openTargets >= 2 && o.sweepTime > 0 && (o.totalTime <= 0 || o.totalTime >= 5));
    for (const auto &s : specs) wanted = wanted || s.name == "sync";
    if (!wanted) return nullptr;
    auto time = std::chrono::steady_clock::now ();
    MatrixCol<T> sumMatrix = MatrixCol<T>::sumProd (-1, pn.getFlowPT (), 1, pn.getFlowTP ());
    // semiflows: non-negative by definition, so every one is a component
    auto [mat, perms] = InvariantMiddle<T>::computePInvariants (sumMatrix, true, std::min<long> (o.timeout, 60),
                                                                o.heuristic (false));
    auto comps = std::make_unique<petri::walk::Components<T>> (wnet, mat);
    std::cout << "Flows: " << mat.getColumnCount () << " P-semiflows in " << millisSince (time) << " ms. ";
    comps->printStats (std::cout);
    return comps;
  }

/**
 * The sweep's choice when the options leave it to us: quests when the net has
 * processes with barriers between them and the targets are conjunctions of
 * place atoms (a fireable is one), rarity otherwise.
 */
template<typename T>
  std::string chooseSweep (const Options &o, const petri::walk::Components<T> *components,
                           const petri::walk::TargetSet<T> &targets)
  {
    if (o.sweepChoice != "auto") return o.sweepChoice;
    if (!components || components->size () == 0 || components->barrierCount () == 0) return "rare";
    for (uint32_t i = 0; i < targets.size (); ++i) {
      const auto &tg = targets.target (i);
      if (!tg.isDeadlock () && !tg.isBound () && tg.expression ().kind != petri::expr::Expression::Kind::Or) return "sync";
    }
    return "rare";
  }

/** The strategy pool of the options: --strategies, else --strategy, else the default pool. */
inline std::vector<petri::walk::StrategySpec> strategyPool (const Options &o)
{
  std::string list = o.strategies;
  if (list.empty ()) {
    // The heuristics saturate: firing the chosen transition while it stays
    // enabled is what crosses a net whose goal sits behind a long repetition
    // of one transition, and it costs nothing where it cannot apply, since
    // maxFirings returns 1 both for a transition that depletes nothing and on
    // a one-safe net. One plain random walk keeps the unsaturated corner.
    list = (o.strategy != "random" || o.threads == 1)
        ? o.strategy : "random,bestfirst+sat,structural+sat,relaxed+sat";
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
          << (st.steps ? st.targetChecks / st.steps : 0) << " checks/step, "
          << st.distinctFired << " distinct transitions fired, " << st.rareEvents << " rare events"
          << (st.policyResets ? ", " + std::to_string (st.policyResets) + " runs ended by the policy" : "")
          << (st.saturations ? ", " + std::to_string (st.saturations) + " saturated" : "") << ")";
      if (rep.strategy != "random") std::cout << ", best heuristic " << rep.minHeuristic;
      if (!rep.notes.empty ()) std::cout << ", " << rep.notes;
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
                const std::vector<petri::walk::StrategySpec> &specs, petri::walk::Knowledge *knowledge,
                const petri::walk::RestartPolicy *policy, const petri::walk::Components<T> *components = nullptr)
  {
    std::unique_ptr<petri::walk::SharedPool<T>> pool;
    if (o.share > 0) pool = std::make_unique<petri::walk::SharedPool<T>> (o.share, o.shareProb);
    std::function<void (const petri::walk::Claim<T>&)> onClaim = [&] (const petri::walk::Claim<T> &c) {
      printFormula (o, targets, c);
    };
    petri::walk::PortfolioResult<T> res = petri::walk::runPortfolio (wnet, targets, focus, specs, o.threads, budget,
                                                                     o.seed, pool.get (), o.debugSteps, onClaim,
                                                                     o.partition, knowledge, policy, components);
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
 * more time would change nothing -- unless --escalate, which raises the step
 * budget tenfold instead and keeps walking while the total budget lasts.
 * A BOUND line is printed each time a walk raised the best value of a bound
 * target, and every bound target ends with one, so a reader that stops
 * listening early still holds the best value known at that time.
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
    std::unique_ptr<petri::walk::Components<T>> components = buildComponents (o, pn, wnet, specs, targets.openCount ());
    std::string sweepChoice = chooseSweep (o, components.get (), targets);
    std::vector<petri::walk::StrategySpec> sweepSpecs = petri::walk::parseStrategySpecs (sweepChoice, 0, 0, o.sample);
    // one memory for the whole file: the sweep and every round add to it
    petri::walk::Knowledge knowledge (wnet.transitionCount ());
    petri::walk::AnyOf sweepPolicy, roundPolicy;
    sweepPolicy.add (std::make_unique<petri::walk::StepBudget> (o.budget.runLength));
    sweepPolicy.add (std::make_unique<petri::walk::WallTime> (o.runTime));
    sweepPolicy.add (std::make_unique<petri::walk::NoveltyStall> (o.noveltyStall));
    roundPolicy.add (std::make_unique<petri::walk::StepBudget> (o.budget.runLength));
    roundPolicy.add (std::make_unique<petri::walk::WallTime> (o.runTime));

    auto walkStart = std::chrono::steady_clock::now ();
    auto elapsedMs = [&] () { return millisSince (walkStart); };
    const long totalMs = o.totalTime * 1000;

    // best value already printed per bound target; a walk that raises one
    // gets it printed at once
    std::vector<long long> printedBound (targets.size (), std::numeric_limits<long long>::min ());
    auto printBounds = [&] (bool all) {
      for (uint32_t i = 0; i < targets.size (); ++i) {
        if (!targets.target (i).isBound ()) continue;
        long long best = targets.bestValue (i);
        if (all || best > printedBound[i]) {
          std::cout << "BOUND " << targets.name (i) << " " << best << std::endl;
          printedBound[i] = best;
        }
      }
    };

    if (o.sweepTime > 0 && targets.openCount () >= 2) {
      long ms = o.sweepTime * 1000;
      if (totalMs > 0) ms = std::min (ms, totalMs);
      std::cout << "Sweep: " << targets.openCount () << " open properties, " << ms << " ms of " << sweepChoice
          << " walks" << (o.sweepChoice == "auto" ? " (chosen)" : "") << ", restarts on " << sweepPolicy.describe () << "."
          << std::endl;
      budget.timeoutMillis = static_cast<uint64_t> (ms);
      runWalk (o, wnet, targets, NO_FOCUS, budget, sweepSpecs, &knowledge, &sweepPolicy, components.get ());
      if (!o.quiet)
        std::cout << "Sweep done: " << knowledge.distinctFired () << " distinct transitions fired by all threads." << std::endl;
      printBounds (false);
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
                                                       targets.target (k).hasHint () ? hintSpecs : specs,
                                                       &knowledge, &roundPolicy, components.get ());
        printBounds (false);
        if (!stepBound (res, budget)) allStepBound = false;
        if (totalMs > 0 && elapsedMs () >= totalMs) break;
      }
      if (totalMs <= 0) break;
      if (allStepBound && targets.openCount () == openBefore && boundValues (targets) == boundsBefore) {
        // Every walk stopped counting steps, not time, and none of them found
        // anything: more time alone would indeed change nothing. More steps
        // might, and under --escalate that is what the remaining budget buys.
        if (!o.escalate || budget.maxSteps == 0) {
          std::cout << "Round " << round << " solved nothing and every walk ended on the step budget: stopping."
              << std::endl;
          break;
        }
        budget.maxSteps *= 10;
        std::cout << "Round " << round << " solved nothing and every walk ended on the step budget: raising it to "
            << budget.maxSteps << " steps." << std::endl;
      }
      perProperty *= 10;
    }
    printBounds (true);
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
    // --strategies, else --strategy, else a pool over the two axes that
    // matter here. Saturation, firing the chosen transition while it stays
    // enabled, is what drains the marking of a net whose dead markings sit
    // behind long repetitions of one transition; the greedy descent of the
    // enabled count is what steers a net where dead markings are rare. Either
    // can be the one that works: on Angiogenesis-PT-20 no unsaturated walk
    // finds the deadlock in ten seconds and every saturated one finds it in
    // milliseconds, while on Philosophers plain random wins the race.
    std::string list = o.strategies;
    if (list.empty ()) {
      list = (o.strategy != "random" || o.threads == 1)
          ? o.strategy : "random,random+sat,deadlock+sat:10:2000,deadlock:30:500";
    }
    std::vector<petri::walk::StrategySpec> specs = petri::walk::parseStrategySpecs (list, o.epsilon, o.stall, o.sample);
    std::vector<petri::walk::Target<T>> tg { petri::walk::Target<T>::deadlockTarget () };
    // A Parikh vector from the caller's state equation describes the firing
    // counts that reach a dead marking, which is the whole answer on a net
    // whose deadlock is structured rather than stumbled into. There is one
    // target here and the caller names it after its own property, so the hint
    // is taken by count, not by name.
    bool hinted = false;
    if (!o.hintsFile.empty ()) {
      auto hints = petri::sexpr::loadHints (o.hintsFile, pn);
      auto it = hints.find ("ReachabilityDeadlock");
      if (it == hints.end () && hints.size () == 1) it = hints.begin ();
      if (it != hints.end ()) {
        tg[0].setHint (it->second);
        hinted = true;
      }
    }
    if (hinted && o.strategies.empty () && o.strategy == "random" && o.threads > 1) {
      specs = petri::walk::parseStrategySpecs ("parikh,deadlock+sat:10:2000,random+sat,random",
                                               o.epsilon, o.stall, o.sample);
    }
    petri::walk::TargetSet<T> targets (pn.getPlaceCount (), std::move (tg), { "ReachabilityDeadlock" }, { "TRUE" });
    petri::walk::Knowledge knowledge (wnet.transitionCount ());
    petri::walk::AnyOf policy;
    policy.add (std::make_unique<petri::walk::StepBudget> (o.budget.runLength));
    policy.add (std::make_unique<petri::walk::WallTime> (o.runTime));
    runWalk (o, wnet, targets, 0, budget, specs, &knowledge, &policy, buildComponents (o, pn, wnet, specs, 0).get ());
  }

} // namespace petri::cli

#endif /* PETRI_CLI_WALKDRIVER_H_ */
