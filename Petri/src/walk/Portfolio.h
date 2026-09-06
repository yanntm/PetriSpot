/*
 * Portfolio.h
 *
 * Several walkers in parallel threads on one TargetSet, each with its own
 * strategy instance, seed, NoveltyTracker and thread-local state, sharing a
 * Knowledge and a RestartPolicy built by the caller. A sweep over at least
 * partitionMin targets is split between the threads, each checking the ids
 * congruent to its rank; below that every thread checks every target. All aim at the same focus
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
#include <limits>
#include <functional>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "walk/BestFirstStrategy.h"
#include "walk/DeadlockStrategy.h"
#include "walk/GoalDistance.h"
#include "walk/ParikhStrategy.h"
#include "walk/RelaxedPlanStrategy.h"
#include "walk/SharedPool.h"
#include "walk/Strategy.h"
#include "walk/StructuralDistance.h"
#include "walk/Components.h"
#include "walk/ComponentStrategy.h"
#include "walk/Knowledge.h"
#include "walk/NoveltyTracker.h"
#include "walk/QuestSweep.h"
#include "walk/RareStrategy.h"
#include "walk/RestartPolicy.h"
#include "walk/TargetSet.h"
#include "walk/Walker.h"

namespace petri::walk
{

struct StrategySpec
{
  std::string name;      // random | rare | bestfirst | structural | relaxed | parikh | sync
  bool saturate = false; // "+sat": fire the chosen transition as many times as possible
  unsigned epsilon = 10; // percent of random moves (heuristic strategies)
  uint64_t stall = 0;    // restart after this many steps without improvement
  size_t sample = 0;     // best-first: candidates scored per step (0: all)

  std::string label () const
  {
    std::string l = name + (saturate ? "+sat" : "");
    if (name == "random" || name == "rare" || name == "parikh") return l;
    return l + ":" + std::to_string (epsilon) + ":" + std::to_string (stall);
  }
};

/** Strip a "+sat" suffix from a strategy name; true iff it was there. */
inline bool stripSaturate (std::string &name)
{
  const std::string suffix = "+sat";
  if (name.size () > suffix.size () && name.compare (name.size () - suffix.size (), suffix.size (), suffix) == 0) {
    name.erase (name.size () - suffix.size ());
    return true;
  }
  return false;
}

/** Parse "name[+sat][:epsilon[:stall]],name..." with defaults for omitted fields. */
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
      if (k == 0) {
        spec.name = field;
        spec.saturate = stripSaturate (spec.name);
      }
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
    DeadlockStrategy<T> *deadlock = nullptr;
    ComponentStrategy<T> *sync = nullptr;
    QuestSweep<T> *sweep = nullptr;

    /** Strategy specific counters for the report; empty when there are none. */
    std::string notes () const
    {
      if (sync)
        return "stages " + std::to_string (sync->replans) + ", barriers fired " + std::to_string (sync->barriersFired)
            + ", depth " + std::to_string (sync->maxDepth) + ", stranded " + std::to_string (sync->hopeless)
            + ", popped " + std::to_string (sync->fallbacks) + ", refusals " + std::to_string (sync->refusals);
      if (sweep)
        return "quests " + std::to_string (sweep->retargets) + ", own target claimed " + std::to_string (sweep->claimedOwn)
            + ", bound steps climbed " + std::to_string (sweep->stepsClimbed) + ", abandoned " + std::to_string (sweep->abandoned)
            + " (hopeless " + std::to_string (sweep->hopelessQuests) + ", unstageable " + std::to_string (sweep->unstageableQuests)
            + "), hopeless markings " + std::to_string (sweep->hopelessMarkings) + ", abandon restarts "
            + std::to_string (sweep->abandonRestarts) + ", filler runs " + std::to_string (sweep->fillerRuns);
      return "";
    }

    /** Best heuristic value seen (instrumentation), or 0 for random. */
    uint64_t minHeuristic () const
    {
      if (sync) return sync->minDistance == std::numeric_limits<uint64_t>::max () ? 0 : sync->minDistance;
      if (bestFirst) return bestFirst->minDistance;
      if (relaxed) return relaxed->minHeuristic;
      if (deadlock) return deadlock->minEnabled == std::numeric_limits<uint64_t>::max () ? 0 : deadlock->minEnabled;
      return 0;
    }
  };

/**
 * The strategy of spec aimed at focus; without a focus (a sweep) the rarity
 * and age choice, or uniform random when the spec says so. A
 * deadlock target has no goal atom either: every heuristic strategy becomes
 * the greedy descent of the enabled count (DeadlockStrategy). A bound without
 * a known limit becomes best-first on the value of the form (BoundDistance).
 */
template<typename T>
  StrategyBundle<T> makeStrategy (const StrategySpec &spec, const WalkNet<T> &net, const Target<T> *focus,
                                  const Components<T> *components = nullptr, const TargetSet<T> *targets = nullptr,
                                  unsigned rank = 0)
  {
    StrategyBundle<T> b;
    b.spec = spec;
    const std::string &n = spec.name;
    if (n == "sync" && !focus && components && targets) {
      // the sweep as quests: the strategy picks its targets itself
      auto s = std::make_unique<QuestSweep<T>> (net, *components, *targets, rank, spec.epsilon, spec.sample, spec.stall);
      b.sweep = s.get ();
      b.strategy = std::move (s);
      return b;
    }
    if (n == "sync" && focus && !focus->isDeadlock () && !focus->isBound () && components) {
      // one process per atom, driven to its place and held there
      auto s = std::make_unique<ComponentStrategy<T>> (net, *components, focus->expression (), spec.epsilon,
                                                        spec.sample, spec.stall);
      if (s->supported ()) {
        b.sync = s.get ();
        b.strategy = std::move (s);
        return b;
      }
      // not a conjunction of place >= k atoms: fall back to the marking distance
      StrategySpec fallback = spec;
      fallback.name = "bestfirst";
      return makeStrategy (fallback, net, focus, components);
    }
    if (n == "sync") {
      StrategySpec fallback = spec;
      fallback.name = focus ? "bestfirst" : "rare";
      return makeStrategy (fallback, net, focus, components);
    }
    if (n == "parikh") {
      if (focus && focus->hasHint ()) {
        b.strategy = std::make_unique<ParikhStrategy<T>> (net, focus->getHint ());
      } else {
        b.strategy = std::make_unique<RandomStrategy<T>> ();
        b.spec.name = "random";
      }
    } else if (focus && focus->isDeadlock ()) {
      // A deadlock target carries no expression, so the distance strategies
      // have nothing to descend; the enabled count is its distance.
      if (n == "random") {
        b.strategy = std::make_unique<RandomStrategy<T>> ();
      } else {
        auto s = std::make_unique<DeadlockStrategy<T>> (net, spec.epsilon, spec.sample, spec.stall);
        b.deadlock = s.get ();
        b.strategy = std::move (s);
        b.spec.name = "deadlock";
      }
    } else if (!focus) {
      if (n == "random") {
        b.strategy = std::make_unique<RandomStrategy<T>> ();
      } else {
        b.strategy = std::make_unique<RareStrategy<T>> (spec.sample);
        b.spec.name = "rare";
      }
    } else if (n == "random") {
      b.strategy = std::make_unique<RandomStrategy<T>> ();
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
  std::string notes;
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
/**
 * Default of partitionMin: 0, sweeps are never split. On ResIsolation-PT-N10P4
 * (95 000 fireable targets) splitting between 4 threads divided the checks per
 * step by four and doubled the step rate, and still claimed 30 % fewer targets
 * in 20 s: a thread walks past the finds another owns. The knob stays for
 * nets where the checks dominate the step.
 */
constexpr size_t PARTITION_MIN = 0;

template<typename T>
  PortfolioResult<T> runPortfolio (const WalkNet<T> &net, TargetSet<T> &targets, uint32_t focus,
                                   const std::vector<StrategySpec> &specs, unsigned threads,
                                   WalkBudget budget, uint64_t seed, SharedPool<T> *pool,
                                   uint64_t debugSteps,
                                   const std::function<void (const Claim<T>&)> &onClaim = {},
                                   size_t partitionMin = PARTITION_MIN, Knowledge *knowledge = nullptr,
                                   const RestartPolicy *restartPolicy = nullptr,
                                   const Components<T> *components = nullptr,
                                   const std::vector<const RestartPolicy*> *policies = nullptr)
  {
    // a policy per strategy of the pool when given (aligned with specs), else the one policy for every thread
    PortfolioResult<T> out;
    out.reports.resize (threads);
    std::atomic<bool> stop (false);
    std::mutex mutex;
    budget.stop = &stop;
    const Target<T> *focusTarget = focus == NO_FOCUS ? nullptr : &targets.target (focus);
    const bool partition = focus == NO_FOCUS && threads > 1 && partitionMin > 0 && targets.size () >= partitionMin;

    auto body = [&] (unsigned i) {
      StrategyBundle<T> bundle = makeStrategy (specs[i % specs.size ()], net, focusTarget, components, &targets, i);
      if (bundle.relaxed && i == 0) bundle.relaxed->debugSteps = debugSteps;
      if (bundle.sweep && i == 0) bundle.sweep->debugSteps = debugSteps;
      if (bundle.sync && i == 0 && debugSteps > 0) {
        bundle.sync->debugSteps = debugSteps;
        bundle.sync->describe (std::cerr, Marking<T> (net.initialMarking ()));
      }
      std::string label = bundle.spec.label ();
      std::vector<uint32_t> subset;
      if (partition) {
        for (uint32_t id = i; id < targets.size (); id += threads) subset.push_back (id);
      }
      Walker<T> walker (net, targets, focus, *bundle.strategy, seed + 7919u * i, partition ? &subset : nullptr);
      NoveltyTracker<T> tracker (net.transitionCount (), knowledge);
      walker.setTracker (&tracker);
      walker.setKnowledge (knowledge);
      walker.setRestartPolicy (policies && !policies->empty () ? (*policies)[i % policies->size ()] : restartPolicy);
      walker.setPool (pool);
      walker.setSaturate (bundle.spec.saturate);
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
      rep.notes = bundle.notes ();
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
