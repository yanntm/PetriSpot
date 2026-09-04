/*
 * Walker.h
 *
 * One explicit walk toward a Target: restart loop, step budget, wall clock,
 * witness trace. See algorithm.md.
 */
#ifndef PETRI_WALK_WALKER_H_
#define PETRI_WALK_WALKER_H_

#include <atomic>
#include <chrono>
#include <cstdint>
#include <random>
#include <vector>

#include "walk/EnabledSet.h"
#include "walk/Marking.h"
#include "walk/SharedPool.h"
#include "walk/Strategy.h"
#include "walk/Target.h"
#include "walk/WalkNet.h"

namespace petri::walk
{

struct WalkStats
{
  uint64_t steps = 0;
  uint64_t resets = 0;
  uint64_t deadEnds = 0;
  uint64_t stalls = 0;    // restarts requested by the strategy
  uint64_t poolRestarts = 0; // restarts drawn from the shared pool
  uint64_t millis = 0;
  uint64_t arcVisits = 0; // consumer arcs visited by the enabled-set update
  uint64_t flips = 0;     // enabled-status changes
};

struct WalkResult
{
  bool found = false;
  bool hasTrace = false;       // trace was recorded (WalkBudget::recordTrace)
  std::vector<uint32_t> trace; // transitions fired from the initial marking
  WalkStats stats;
};

struct WalkBudget
{
  uint64_t maxSteps = 0;       // 0: unlimited
  uint64_t runLength = 1000000; // steps between resets
  uint64_t timeoutMillis = 0;  // 0: unlimited
  bool recordTrace = false;    // keep the fired transitions of the current run
  const std::atomic<bool> *stop = nullptr; // external stop request, polled
};

template<typename T>
  class Walker
  {
    const WalkNet<T> &net;
    const Target<T> &target;
    Strategy<T> &strategy;
    std::mt19937_64 rng;

    Marking<T> initialMarking;
    EnabledSet<T> initialEnabled;
    Marking<T> marking;
    EnabledSet<T> enabled;
    std::vector<uint32_t> trace;
    bool recordTrace = false;

    SharedPool<T> *pool = nullptr;
    typename SharedPool<T>::Entry poolEntry; // origin of the current run when drawn
    bool fromPool = false;
    SparseArray<T> scratchMarking;

    /** Offer the best state of the finished run to the pool, report the origin. */
    void publishRun ()
    {
      if (!pool) return;
      uint64_t h = 0;
      bool has = strategy.bestOfRun (scratchMarking, h);
      if (fromPool) pool->report (poolEntry.id, has && h < poolEntry.heuristic);
      if (has && (!fromPool || h < poolEntry.heuristic)) {
        // the trace to the published state is not tracked; share it only
        // when traces are off (a pooled restart then yields no witness trace)
        pool->publish (scratchMarking, h, std::vector<uint32_t> ());
      }
    }

    void reset (WalkStats *st = nullptr)
    {
      publishRun ();
      fromPool = false;
      if (pool && pool->draw (rng, poolEntry)) {
        fromPool = true;
        marking = Marking<T> (poolEntry.marking);
        enabled.initialize (marking);
        if (st) ++st->poolRestarts;
      } else {
        marking.assign (initialMarking);
        enabled.assign (initialEnabled);
      }
      trace.clear ();
      strategy.onReset ();
    }

    void fire (uint32_t t)
    {
      marking.apply (net.effect (t), [this] (size_t p, T oldv, T newv) {
        enabled.onPlaceChanged (p, oldv, newv);
      });
      if (recordTrace) trace.push_back (t);
    }

  public:
    Walker (const WalkNet<T> &n, const Target<T> &tg, Strategy<T> &st,
            uint64_t seed)
        : net (n), target (tg), strategy (st), rng (seed),
          initialMarking (n.initialMarking ()), initialEnabled (n),
          marking (n.initialMarking ()), enabled (n)
    {
      initialEnabled.initialize (initialMarking);
      enabled.assign (initialEnabled);
    }

    WalkResult run (const WalkBudget &budget)
    {
      using clock = std::chrono::steady_clock;
      auto start = clock::now ();
      WalkResult result;
      recordTrace = budget.recordTrace;
      result.hasTrace = recordTrace;
      WalkStats &st = result.stats;
      WalkContext<T> ctx { net, marking, enabled, rng };
      uint64_t runSteps = 0;
      fromPool = false;
      marking.assign (initialMarking);
      enabled.assign (initialEnabled);
      trace.clear ();
      strategy.onReset ();
      for (;;) {
        if (target.reached (marking, enabled)) {
          result.found = true;
          if (recordTrace && !fromPool) result.trace = trace;
          else result.hasTrace = false;
          break;
        }
        if (enabled.empty () || runSteps >= budget.runLength) {
          if (enabled.empty ()) ++st.deadEnds;
          ++st.resets;
          runSteps = 0;
          reset (&st);
          continue;
        }
        if ((st.steps & 1023) == 0) {
          if (budget.stop && budget.stop->load (std::memory_order_relaxed)) break;
          if (budget.maxSteps && st.steps >= budget.maxSteps) break;
          if (budget.timeoutMillis) {
            auto ms = std::chrono::duration_cast<std::chrono::milliseconds> (clock::now () - start).count ();
            if (static_cast<uint64_t> (ms) >= budget.timeoutMillis) break;
          }
        }
        uint32_t t = strategy.choose (ctx);
        if (t == RESTART) {
          ++st.stalls;
          ++st.resets;
          runSteps = 0;
          reset (&st);
          continue;
        }
        fire (t);
        ++st.steps;
        ++runSteps;
      }
      st.millis = static_cast<uint64_t> (
          std::chrono::duration_cast<std::chrono::milliseconds> (clock::now () - start).count ());
      st.arcVisits = enabled.arcVisits;
      st.flips = enabled.flips;
      return result;
    }

    /** Replay a trace from the initial marking; true iff the target holds at its end. */
    bool verify (const std::vector<uint32_t> &witness)
    {
      recordTrace = false;
      SharedPool<T> *saved = pool;
      pool = nullptr;
      reset ();
      pool = saved;
      for (uint32_t t : witness) {
        if (!enabled.isEnabled (t)) return false;
        fire (t);
      }
      return target.reached (marking, enabled);
    }

    const Marking<T>& currentMarking () const
    {
      return marking;
    }
    void setPool (SharedPool<T> *p)
    {
      pool = p;
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_WALKER_H_ */
