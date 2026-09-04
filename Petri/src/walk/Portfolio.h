/*
 * Portfolio.h
 *
 * Several walkers in parallel threads on one TargetSet, each with its own
 * strategy instance, seed and thread-local state. All aim at the same focus
 * (or none) and claim any target they reach; claims are published in order
 * through a callback, the focus claim stops the others, recorded traces are
 * verified before being reported. Strategies are named by specs
 * ("relaxed:0:300" = name:epsilon:stall) assigned round-robin to threads. An
 * optional SharedPool lets restarts draw promising states from other threads.
 */
#ifndef PETRI_WALK_PORTFOLIO_H_
#define PETRI_WALK_PORTFOLIO_H_

#include <atomic>
#include <cstdint>
#include <functional>
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
#include "walk/TargetSet.h"
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

/**
 * The strategy of spec aimed at focus; random when there is no focus or it is
 * a deadlock. A bound without a known limit has no goal atom: every heuristic
 * strategy becomes best-first on the value of the form (BoundDistance).
 */
template<typename T>
  StrategyBundle<T> makeStrategy (const StrategySpec &spec, const WalkNet<T> &net, const Target<T> *focus)
  {
    StrategyBundle<T> b;
    b.spec = spec;
    const std::string &n = spec.name;
    if (!focus || focus->isDeadlock () || n == "random") {
      b.strategy = std::make_unique<RandomStrategy<T>> ();
      b.spec.name = "random";
    } else if (focus->isBound () && !focus->hasLimit ()) {
      b.distance = std::make_unique<BoundDistance<T>> (focus->boundForm ());
      auto s = std::make_unique<BestFirstStrategy<T>> (*b.distance, spec.epsilon, spec.sample, spec.stall);
      b.bestFirst = s.get ();
      b.strategy = std::move (s);
      b.spec.name = "bestfirst";
    } else if (n == "bestfirst" || n == "structural") {
      if (n == "bestfirst") b.distance = std::make_unique<MarkingDistance<T>> (focus->expression ());
      else b.distance = std::make_unique<StructuralDistance<T>> (focus->expression (), net);
      auto s = std::make_unique<BestFirstStrategy<T>> (*b.distance, spec.epsilon, spec.sample, spec.stall);
      b.bestFirst = s.get ();
      b.strategy = std::move (s);
    } else if (n == "relaxed") {
      auto s = std::make_unique<RelaxedPlanStrategy<T>> (net, focus->expression (), spec.epsilon, spec.stall);
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
  bool found = false;   // claimed the focus
  uint64_t claims = 0;
};

/** A target won by a thread: where, by whom, and how it was reached. */
template<typename T>
  struct Claim
  {
    uint32_t target = 0;
    unsigned thread = 0;
    std::string strategy;       // label of the claiming thread's strategy
    SparseArray<T> marking;     // marking reached
    bool hasTrace = false;      // trace recorded and verified
    std::vector<uint32_t> trace;
  };

template<typename T>
  struct PortfolioResult
  {
    bool found = false;         // the focus was claimed during this run
    size_t winner = 0;          // thread that claimed the focus
    std::string winnerStrategy;
    std::vector<Claim<T>> claims; // every target claimed, in claim order
    std::vector<ThreadReport> reports;
  };

/**
 * Run threads walkers on targets toward focus (NO_FOCUS: sweep). onClaim, if
 * given, is called under a lock as each claim is published; traces are only
 * verified after the threads join (Claim::hasTrace).
 */
template<typename T>
  PortfolioResult<T> runPortfolio (const WalkNet<T> &net, TargetSet<T> &targets, uint32_t focus,
                                   const std::vector<StrategySpec> &specs, unsigned threads,
                                   WalkBudget budget, uint64_t seed, SharedPool<T> *pool,
                                   uint64_t debugSteps,
                                   const std::function<void (const Claim<T>&)> &onClaim = {})
  {
    PortfolioResult<T> out;
    out.reports.resize (threads);
    std::atomic<bool> stop (false);
    std::mutex mutex;
    budget.stop = &stop;
    const Target<T> *focusTarget = focus == NO_FOCUS ? nullptr : &targets.target (focus);

    auto body = [&] (unsigned i) {
      StrategyBundle<T> bundle = makeStrategy (specs[i % specs.size ()], net, focusTarget);
      if (bundle.relaxed && i == 0) bundle.relaxed->debugSteps = debugSteps;
      std::string label = bundle.spec.label ();
      Walker<T> walker (net, targets, focus, *bundle.strategy, seed + 7919u * i);
      walker.setPool (pool);
      walker.setOnClaim ([&, i, label] (uint32_t id, const Marking<T> &m, const std::vector<uint32_t> *trace) {
        Claim<T> c;
        c.target = id;
        c.thread = i;
        c.strategy = label;
        c.marking = m.sparse ();
        if (trace) {
          c.hasTrace = true;
          c.trace = *trace;
        }
        std::lock_guard<std::mutex> lock (mutex);
        if (id == focus) {
          out.found = true;
          out.winner = i;
          out.winnerStrategy = label;
          stop.store (true);
        }
        out.claims.push_back (std::move (c));
        if (onClaim) onClaim (out.claims.back ());
      });
      WalkResult res = walker.run (budget);
      ThreadReport &rep = out.reports[i];
      rep.strategy = label;
      rep.stats = res.stats;
      rep.minHeuristic = bundle.minHeuristic ();
      rep.found = res.found;
      rep.claims = res.claims;
    };

    if (threads <= 1) {
      body (0);
    } else {
      std::vector<std::thread> pool_;
      for (unsigned i = 0; i < threads; ++i) pool_.emplace_back (body, i);
      for (auto &th : pool_) th.join ();
    }

    bool anyTrace = false;
    for (const auto &c : out.claims) anyTrace = anyTrace || c.hasTrace;
    if (anyTrace) {
      RandomStrategy<T> rs;
      Walker<T> verifier (net, targets, NO_FOCUS, rs, 0);
      for (auto &c : out.claims) {
        if (c.hasTrace && !verifier.verify (c.target, c.trace)) {
          std::cerr << "Internal error: witness trace of thread " << c.thread << " for " << targets.name (c.target)
              << " does not replay to the goal." << std::endl;
          c.hasTrace = false;
        }
      }
    }
    return out;
  }

} // namespace petri::walk

#endif /* PETRI_WALK_PORTFOLIO_H_ */
