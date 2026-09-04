/*
 * WalkDriver.h
 *
 * Reachability from the command line: load and print properties, schedule
 * the open ones in rounds, run the portfolio on one target and print the
 * MCC verdict line.
 */
#ifndef PETRI_CLI_WALKDRIVER_H_
#define PETRI_CLI_WALKDRIVER_H_

#include <algorithm>
#include <chrono>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "cli/Options.h"
#include "core/SparsePetriNet.h"
#include "expr/Property.h"
#include "expr/SexprPrinter.h"
#include "parse/PropertyFile.h"
#include "walk/Portfolio.h"
#include "walk/SharedPool.h"
#include "walk/Target.h"
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

/** Run the portfolio toward target; print the MCC verdict line if reached. */
template<typename T>
  bool runWalk (const Options &o, const petri::walk::WalkNet<T> &wnet, const petri::walk::Target<T> &target,
                const std::string &name, const std::string &verdict, const petri::walk::WalkBudget &budget,
                const std::vector<petri::walk::StrategySpec> &specs)
  {
    std::unique_ptr<petri::walk::SharedPool<T>> pool;
    if (o.share > 0) pool = std::make_unique<petri::walk::SharedPool<T>> (o.share, o.shareProb);
    petri::walk::PortfolioResult<T> res = petri::walk::runPortfolio (wnet, target, specs, o.threads, budget,
                                                                     o.seed, pool.get (), o.debugSteps);
    for (size_t i = 0; i < res.reports.size (); ++i) {
      const petri::walk::ThreadReport &rep = res.reports[i];
      const petri::walk::WalkStats &st = rep.stats;
      uint64_t ms = st.millis == 0 ? 1 : st.millis;
      if (o.quiet && !rep.found && res.found) continue;
      std::cout << "Thread " << i << " [" << rep.strategy << "] " << (rep.found ? "found a witness" : "stopped")
          << " after " << st.steps << " steps, " << st.resets << " resets (" << st.deadEnds << " dead ends, "
          << st.stalls << " stalls, " << st.poolRestarts << " from pool), " << st.millis << " ms ("
          << (st.steps / ms) << " steps/ms; " << (st.steps ? st.arcVisits / st.steps : 0) << " arc visits/step)";
      if (rep.strategy != "random") std::cout << ", best heuristic " << rep.minHeuristic;
      std::cout << "." << std::endl;
    }
    if (pool) {
      std::cout << "Shared pool: " << pool->publishedCount () << " published, " << pool->drawnCount ()
          << " drawn, " << pool->evictedCount () << " evicted, " << pool->size () << " held." << std::endl;
    }
    if (!res.found) {
      std::cout << "No witness found for " << name << "." << std::endl;
      return false;
    }
    std::cout << "FORMULA " << name << " " << verdict << " TECHNIQUES EXPLICIT "
        << (res.winnerStrategy == "random" ? "RANDOM_WALK" : "HEURISTIC_WALK")
        << (o.threads > 1 ? " PARALLEL_PROCESSING" : "") << std::endl;
    if (res.result.hasTrace) {
      const auto &tnames = wnet.getNet ().getTnames ();
      std::cout << "Witness (" << res.result.trace.size () << " transitions):";
      for (uint32_t t : res.result.trace) std::cout << " " << tnames[t];
      std::cout << std::endl;
    } else if (!o.quiet) {
      std::cout << "Witness marking ";
      printMarking (std::cout, res.witness, wnet.getNet ().getPnames ());
      std::cout << std::endl;
    }
    return true;
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

/**
 * Walk every property of the file: trivial ones are answered at once, the
 * others are scheduled in rounds of growing per-property budget when
 * --totalTime is set, or each get the -t timeout otherwise.
 */
template<typename T>
  void runProperties (const Options &o, const SparsePetriNet<T> &pn)
  {
    using petri::expr::PropertyKind;
    std::vector<petri::expr::Property> props = loadProperties (o, pn);
    if (o.printProps) {
      printProperties (props, pn, o.printPropsFormat);
      return;
    }
    petri::walk::WalkBudget budget = o.budget;
    budget.recordTrace = o.trace;
    petri::walk::WalkNet<T> wnet (pn);
    std::vector<petri::walk::StrategySpec> specs = strategyPool (o);
    std::vector<size_t> open;
    for (size_t k = 0; k < props.size (); ++k) {
      const auto &prop = props[k];
      if (prop.kind == PropertyKind::Unsupported) {
        std::cout << "Skipping " << prop.name << " : " << prop.comment << std::endl;
        continue;
      }
      if (prop.kind != PropertyKind::Deadlock && prop.goal ().kind == petri::expr::Expression::Kind::False) {
        std::cout << "FORMULA " << prop.name << " " << (prop.kind == PropertyKind::Invariant ? "TRUE" : "FALSE")
            << " TECHNIQUES TOPOLOGICAL TRIVIAL" << std::endl;
        continue;
      }
      open.push_back (k);
    }
    auto walkStart = std::chrono::steady_clock::now ();
    auto elapsedSeconds = [&] () {
      return std::chrono::duration_cast<std::chrono::seconds> (std::chrono::steady_clock::now () - walkStart).count ();
    };
    long perProperty = o.totalTime > 0 ? o.roundTime : o.timeout;
    for (int round = 1; !open.empty (); ++round) {
      if (o.totalTime > 0) {
        long left = o.totalTime - elapsedSeconds ();
        if (left < static_cast<long> (open.size ())) break; // not even a second each
        // last round when the geometric budget would exceed what is left
        if (perProperty * static_cast<long> (open.size ()) >= left)
          perProperty = std::max (1L, left / static_cast<long> (open.size ()));
        std::cout << "Round " << round << ": " << open.size () << " open properties, " << perProperty
            << " s each, " << left << " s left." << std::endl;
      }
      budget.timeoutMillis = static_cast<uint64_t> (perProperty) * 1000;
      std::vector<size_t> still;
      for (size_t i = 0; i < open.size (); ++i) {
        const auto &prop = props[open[i]];
        petri::walk::Target<T> target = prop.kind == PropertyKind::Deadlock
            ? petri::walk::Target<T>::deadlockTarget () : petri::walk::Target<T> (prop.goal ());
        bool solved = runWalk (o, wnet, target, prop.name, prop.verdictIfReached (), budget, specs);
        if (!solved) still.push_back (open[i]);
        if (o.totalTime > 0 && elapsedSeconds () >= o.totalTime) {
          still.insert (still.end (), open.begin () + static_cast<long> (i) + 1, open.end ());
          break;
        }
      }
      open = still;
      if (o.totalTime <= 0) break;
      perProperty *= 10;
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
    petri::walk::Target<T> target = petri::walk::Target<T>::deadlockTarget ();
    runWalk (o, wnet, target, "ReachabilityDeadlock", "TRUE", budget, specs);
  }

} // namespace petri::cli

#endif /* PETRI_CLI_WALKDRIVER_H_ */
