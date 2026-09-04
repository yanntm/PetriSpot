/*
 * Walker.h
 *
 * One explicit walk over a TargetSet with an optional focus: restart loop,
 * step budget, wall clock, incremental target checks, claims, witness trace.
 * See algorithm.md.
 */
#ifndef PETRI_WALK_WALKER_H_
#define PETRI_WALK_WALKER_H_

#include <atomic>
#include <chrono>
#include <cstdint>
#include <functional>
#include <random>
#include <vector>

#include "walk/EnabledSet.h"
#include "walk/Marking.h"
#include "walk/SharedPool.h"
#include "walk/Strategy.h"
#include "walk/TargetSet.h"
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
  uint64_t targetChecks = 0; // goal evaluations
};

struct WalkResult
{
  bool found = false;    // the focus was claimed by this walker
  uint64_t claims = 0;   // targets claimed by this walker, focus included
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
  public:
    /**
     * Called when this walker wins a target: the marking reached and the
     * trace from the initial marking, or nullptr when not recorded or when
     * the run started from a pooled state.
     */
    using OnClaim = std::function<void (uint32_t target, const Marking<T> &marking,
                                        const std::vector<uint32_t> *trace)>;

  private:
    const WalkNet<T> &net;
    TargetSet<T> &targets;
    uint32_t focus;
    Strategy<T> &strategy;
    std::mt19937_64 rng;
    OnClaim onClaim;

    Marking<T> initialMarking;
    EnabledSet<T> initialEnabled;
    Marking<T> marking;
    EnabledSet<T> enabled;
    std::vector<uint32_t> trace;
    bool recordTrace = false;

    std::vector<size_t> touched;    // places changed by the last firing
    std::vector<uint32_t> stamp;    // per target: epoch of its last candidacy
    std::vector<uint32_t> candidates;
    uint32_t epoch = 0;

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
      touched.clear ();
      marking.apply (net.effect (t), [this] (size_t p, T oldv, T newv) {
        enabled.onPlaceChanged (p, oldv, newv);
        touched.push_back (p);
      });
      if (recordTrace) trace.push_back (t);
    }

    /** Evaluate an open target here; claim it if it holds. */
    void check (uint32_t id, WalkResult &result)
    {
      if (targets.isSolved (id)) return;
      ++result.stats.targetChecks;
      if (!targets.target (id).reached (marking, enabled)) return;
      if (!targets.claim (id)) return;
      ++result.claims;
      if (id == focus) result.found = true;
      if (onClaim) onClaim (id, marking, recordTrace && !fromPool ? &trace : nullptr);
    }

    void checkAll (WalkResult &result)
    {
      for (uint32_t id = 0; id < targets.size (); ++id)
        if (!targets.target (id).isDeadlock ()) check (id, result);
    }

    /** The open targets mentioning a place changed by the last firing. */
    void checkTouched (WalkResult &result)
    {
      ++epoch;
      candidates.clear ();
      for (size_t p : touched) {
        for (uint32_t id : targets.targetsOf (p)) {
          if (stamp[id] == epoch) continue;
          stamp[id] = epoch;
          candidates.push_back (id);
        }
      }
      for (uint32_t id : candidates) check (id, result);
    }

    void checkDeadlocks (WalkResult &result)
    {
      for (uint32_t id : targets.deadlocks ()) check (id, result);
    }

    bool done () const
    {
      if (focus != NO_FOCUS && targets.isSolved (focus)) return true;
      return targets.openCount () == 0;
    }

  public:
    Walker (const WalkNet<T> &n, TargetSet<T> &tgs, uint32_t focusTarget, Strategy<T> &st, uint64_t seed)
        : net (n), targets (tgs), focus (focusTarget), strategy (st), rng (seed),
          initialMarking (n.initialMarking ()), initialEnabled (n),
          marking (n.initialMarking ()), enabled (n), stamp (tgs.size (), 0)
    {
      initialEnabled.initialize (initialMarking);
      enabled.assign (initialEnabled);
    }

    void setOnClaim (OnClaim cb)
    {
      onClaim = std::move (cb);
    }
    void setPool (SharedPool<T> *p)
    {
      pool = p;
    }

    WalkResult run (const WalkBudget &budget)
    {
      using clock = std::chrono::steady_clock;
      auto start = clock::now ();
      WalkResult result;
      recordTrace = budget.recordTrace;
      WalkStats &st = result.stats;
      WalkContext<T> ctx { net, marking, enabled, rng };
      uint64_t runSteps = 0;
      fromPool = false;
      marking.assign (initialMarking);
      enabled.assign (initialEnabled);
      trace.clear ();
      strategy.onReset ();
      checkAll (result);
      while (!done ()) {
        if (enabled.empty ()) {
          checkDeadlocks (result);
          if (done ()) break;
          ++st.deadEnds;
        }
        if (enabled.empty () || runSteps >= budget.runLength) {
          ++st.resets;
          runSteps = 0;
          reset (&st);
          checkAll (result);
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
          checkAll (result);
          continue;
        }
        fire (t);
        ++st.steps;
        ++runSteps;
        checkTouched (result);
      }
      st.millis = static_cast<uint64_t> (
          std::chrono::duration_cast<std::chrono::milliseconds> (clock::now () - start).count ());
      st.arcVisits = enabled.arcVisits;
      st.flips = enabled.flips;
      return result;
    }

    /** Replay a trace from the initial marking; true iff target id holds at its end. */
    bool verify (uint32_t id, const std::vector<uint32_t> &witness)
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
      return targets.target (id).reached (marking, enabled);
    }

    const Marking<T>& currentMarking () const
    {
      return marking;
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_WALKER_H_ */
