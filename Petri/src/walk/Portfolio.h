/*
 * Portfolio.h
 *
 * Several walkers in parallel threads, each with its own strategy instance,
 * seed and thread-local state, racing toward the same Target. The first
 * verified witness stops the others. Strategies are named by specs
 * ("relaxed:0:300" = name:epsilon:stall) assigned round-robin to threads. An
 * optional SharedPool lets restarts draw promising states from other threads.
 */
#ifndef PETRI_WALK_PORTFOLIO_H_
#define PETRI_WALK_PORTFOLIO_H_

#include <atomic>
#include <cstdint>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "walk/BestFirstStrategy.h"
#include "walk/GoalDistance.h"
#include "walk/RelaxedPlanStrategy.h"
#include "walk/SharedPool.h"
#include "walk/Strategy.h"
#include "walk/StructuralDistance.h"
#include "walk/Target.h"
#include "walk/Walker.h"

namespace petri::walk
{

struct StrategySpec
{
  std::string name;      // random | bestfirst | structural | relaxed
  unsigned epsilon = 10; // percent of random moves (heuristic strategies)
  uint64_t stall = 0;    // restart after this many steps without improvement
  size_t sample = 0;     // best-first: candidates scored per step (0: all)

  std::string label () const
  {
    if (name == "random") return name;
    return name + ":" + std::to_string (epsilon) + ":" + std::to_string (stall);
  }
};

/** Parse "name[:epsilon[:stall]],name..." with defaults for omitted fields. */
inline std::vector<StrategySpec> parseStrategySpecs (const std::string &list,
                                                     unsigned epsilon, uint64_t stall,
                                                     size_t sample)
{
  std::vector<StrategySpec> specs;
  std::stringstream ss (list);
  std::string item;
  while (std::getline (ss, item, ',')) {
    if (item.empty ()) continue;
    StrategySpec spec;
    spec.epsilon = epsilon;
    spec.stall = stall;
    spec.sample = sample;
    std::stringstream is (item);
    std::string field;
    int k = 0;
    while (std::getline (is, field, ':')) {
      if (k == 0) spec.name = field;
      else if (k == 1 && !field.empty ()) spec.epsilon = static_cast<unsigned> (std::stoul (field));
      else if (k == 2 && !field.empty ()) spec.stall = std::stoull (field);
      ++k;
    }
    specs.push_back (spec);
  }
  return specs;
}

/** A strategy instance and the objects it depends on. */
template<typename T>
  struct StrategyBundle
  {
    StrategySpec spec;
    std::unique_ptr<GoalDistance<T>> distance;
    std::unique_ptr<Strategy<T>> strategy;
    BestFirstStrategy<T> *bestFirst = nullptr;
    RelaxedPlanStrategy<T> *relaxed = nullptr;

    /** Best heuristic value seen (instrumentation), or 0 for random. */
    uint64_t minHeuristic () const
    {
      if (bestFirst) return bestFirst->minDistance;
      if (relaxed) return relaxed->minHeuristic;
      return 0;
    }
  };

template<typename T>
  StrategyBundle<T> makeStrategy (const StrategySpec &spec, const WalkNet<T> &net,
                                  const Target<T> &target)
  {
    StrategyBundle<T> b;
    b.spec = spec;
    const std::string &n = spec.name;
    if (target.isDeadlock () || n == "random") {
      b.strategy = std::make_unique<RandomStrategy<T>> ();
      b.spec.name = "random";
    } else if (n == "bestfirst" || n == "structural") {
      if (n == "bestfirst") b.distance = std::make_unique<MarkingDistance<T>> (target.expression ());
      else b.distance = std::make_unique<StructuralDistance<T>> (target.expression (), net);
      auto s = std::make_unique<BestFirstStrategy<T>> (*b.distance, spec.epsilon, spec.sample, spec.stall);
      b.bestFirst = s.get ();
      b.strategy = std::move (s);
    } else if (n == "relaxed") {
      auto s = std::make_unique<RelaxedPlanStrategy<T>> (net, target.expression (), spec.epsilon, spec.stall);
      b.relaxed = s.get ();
      b.strategy = std::move (s);
    } else {
      throw std::string ("Unknown strategy: " + n);
    }
    return b;
  }

struct ThreadReport
{
  std::string strategy;
  WalkStats stats;
  uint64_t minHeuristic = 0;
  bool found = false;
};

template<typename T>
  struct PortfolioResult
  {
    bool found = false;
    size_t winner = 0;
    std::string winnerStrategy;
    WalkResult result;          // of the winner (trace only if recorded and verified)
    SparseArray<T> witness;     // marking reached by the winner
    std::vector<ThreadReport> reports;
  };

template<typename T>
  PortfolioResult<T> runPortfolio (const WalkNet<T> &net, const Target<T> &target,
                                   const std::vector<StrategySpec> &specs, unsigned threads,
                                   WalkBudget budget, uint64_t seed, SharedPool<T> *pool,
                                   uint64_t debugSteps)
  {
    PortfolioResult<T> out;
    out.reports.resize (threads);
    std::atomic<bool> stop (false);
    std::mutex mutex;
    budget.stop = &stop;

    auto body = [&] (unsigned i) {
      StrategyBundle<T> bundle = makeStrategy (specs[i % specs.size ()], net, target);
      if (bundle.relaxed && i == 0) bundle.relaxed->debugSteps = debugSteps;
      Walker<T> walker (net, target, *bundle.strategy, seed + 7919u * i);
      walker.setPool (pool);
      WalkResult res = walker.run (budget);
      ThreadReport &rep = out.reports[i];
      rep.strategy = bundle.spec.label ();
      rep.stats = res.stats;
      rep.minHeuristic = bundle.minHeuristic ();
      rep.found = res.found;
      if (!res.found) return;
      SparseArray<T> witness = walker.currentMarking ().sparse ();
      if (res.hasTrace && !walker.verify (res.trace)) {
        std::lock_guard<std::mutex> lock (mutex);
        std::cerr << "Internal error: witness trace of thread " << i << " does not replay to the goal." << std::endl;
        return;
      }
      std::lock_guard<std::mutex> lock (mutex);
      if (!out.found) {
        out.found = true;
        out.winner = i;
        out.winnerStrategy = bundle.spec.label ();
        out.result = std::move (res);
        out.witness = std::move (witness);
        stop.store (true);
      }
    };

    if (threads <= 1) {
      body (0);
    } else {
      std::vector<std::thread> pool_;
      for (unsigned i = 0; i < threads; ++i) pool_.emplace_back (body, i);
      for (auto &th : pool_) th.join ();
    }
    return out;
  }

} // namespace petri::walk

#endif /* PETRI_WALK_PORTFOLIO_H_ */
