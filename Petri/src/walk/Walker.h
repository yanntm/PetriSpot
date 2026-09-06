/*
 * Walker.h
 *
 * One explicit walk over a TargetSet with an optional focus: the step loop,
 * the restarts a RestartPolicy asks for, incremental target checks over a
 * thread-local index of the targets it owns (all of them, or a subset of a
 * large set), claims, witness trace. A NoveltyTracker, when injected, keeps
 * the walker's firing memory and feeds the shared Knowledge and the pool.
 * See algorithm.md.
 */
#ifndef PETRI_WALK_WALKER_H_
#define PETRI_WALK_WALKER_H_

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <functional>
#include <limits>
#include <random>
#include <vector>

#include "walk/EnabledSet.h"
#include "walk/Task.h"
#include "walk/Knowledge.h"
#include "walk/Marking.h"
#include "walk/NoveltyTracker.h"
#include "walk/RestartPolicy.h"
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
  uint64_t condemned = 0;    // pooled start states found hopeless on arrival and evicted
  uint64_t millis = 0;
  uint64_t arcVisits = 0; // consumer arcs visited by the enabled-set update
  uint64_t flips = 0;     // enabled-status changes
  uint64_t targetChecks = 0; // goal evaluations
  uint64_t saturations = 0;  // updates that fired a transition more than once
  uint64_t policyResets = 0; // runs ended by the restart policy (the others: dead ends, stalls)
  uint64_t inPlaceResets = 0; // of which the marking was kept, the strategy alone restarting
  uint64_t distinctFired = 0; // transitions fired at least once, a coverage measure of the walk
  uint64_t rareEvents = 0;   // first firings, by anyone, of a transition
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
  uint64_t runLength = 1000000; // steps between resets, when no RestartPolicy is injected
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
    bool saturate = false;          // fire the chosen transition as many times as the marking allows

    std::vector<std::pair<size_t, bool>> touched; // places changed by the last firing, and whether they grew
    std::vector<long long> published; // per bound target: last value sent to the shared maximum
    std::vector<uint32_t> stamp;    // per target: epoch of its last candidacy
    std::vector<uint32_t> candidates;
    uint32_t epoch = 0;

    // The targets this walker checks, and per place the ones an increase (up)
    // or a decrease (down) of the place can make hold. Solved ids are dropped
    // from these lists as they are met, so the lists shrink over the run.
    std::vector<uint32_t> own;
    std::vector<std::vector<uint32_t>> up, down;
    std::vector<typename TargetSet<T>::Mention> mentionBuf;

    // Injected memory and policy; both optional.
    NoveltyTracker<T> *tracker = nullptr;
    const Knowledge *knowledge = nullptr;
    const RestartPolicy *restartPolicy = nullptr;

    // The resumable loop's state between two slices (begin / runSlice / finish).
    WalkBudget budget;
    WalkResult result;
    StepBudget fallback { 1000000 };
    const RestartPolicy *policy = nullptr;
    const AnyOf *composite = nullptr;
    uint64_t runSteps = 0, iterations = 0, activeMs = 0, runStartMs = 0;
    bool over = false;

    /** Index the given targets; deadlock targets are checked apart and skipped here. */
    void buildIndex (const std::vector<uint32_t> &ids)
    {
      own.clear ();
      for (auto &l : up) l.clear ();
      for (auto &l : down) l.clear ();
      for (uint32_t id : ids) {
        if (targets.target (id).isDeadlock () || targets.isSolved (id)) continue;
        own.push_back (id);
        targets.mentions (id, mentionBuf);
        for (const auto &m : mentionBuf) {
          if (m.direction >= 0) up[m.place].push_back (id);
          if (m.direction <= 0) down[m.place].push_back (id);
        }
      }
    }

    /** A walker whose own targets are all solved takes over every open one. */
    void refill ()
    {
      if (own.empty () && targets.openCount () > targets.deadlocks ().size ()) buildIndex (targets.openTargets ());
    }

    /** Drop list[i] as solved; the caller does not advance. */
    static void dropAt (std::vector<uint32_t> &list, size_t i)
    {
      list[i] = list.back ();
      list.pop_back ();
    }

    SharedPool<T> *pool = nullptr;
    typename SharedPool<T>::Entry poolEntry; // origin of the current run when drawn
    bool fromPool = false;
    SparseArray<T> scratchMarking;

    /** Offer the best state of the finished run and the last rare event to the pool, report the origin. */
    void publishRun ()
    {
      if (!pool) return;
      if (tracker && tracker->takeRareMarking (scratchMarking)) {
        // a marking nobody had reached the way to: worth a restart, whatever its heuristic
        pool->publish (scratchMarking, 0, std::vector<uint32_t> ());
      }
      uint64_t h = 0;
      bool has = strategy.bestOfRun (scratchMarking, h);
      if (fromPool) pool->report (poolEntry.id, has && h < poolEntry.heuristic);
      if (has && (!fromPool || h < poolEntry.heuristic)) {
        // the trace to the published state is not tracked; share it only
        // when traces are off (a pooled restart then yields no witness trace)
        pool->publish (scratchMarking, h, std::vector<uint32_t> ());
      }
    }

    /** The strategy forgets, the walk goes on from the current marking. */
    void resetInPlace ()
    {
      publishRun ();
      if (tracker) {
        tracker->merge ();
        tracker->resetNovelty ();
      }
      strategy.onReset ();
    }

    void reset (WalkStats *st = nullptr)
    {
      publishRun ();
      if (tracker) {
        tracker->merge ();
        tracker->resetNovelty ();
      }
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

    /**
     * How many times t can be fired in a row from the current marking: it stays
     * enabled as long as every place it depletes keeps its input weight;
     * 1 when t depletes nothing (it would never stop).
     */
    uint64_t maxFirings (uint32_t t) const
    {
      const SparseArray<T> &eff = net.effect (t);
      const SparseArray<T> &pre = net.pre (t);
      uint64_t k = std::numeric_limits<uint64_t>::max ();
      for (size_t i = 0; i < eff.size (); ++i) {
        T d = eff.valueAt (i);
        if (d >= 0) continue;
        size_t p = eff.keyAt (i);
        T have = marking.get (p) - pre.get (p); // tokens beyond the first firing
        k = std::min (k, static_cast<uint64_t> (have / -d) + 1);
      }
      return k == std::numeric_limits<uint64_t>::max () ? 1 : k;
    }

    void fire (uint32_t t, uint64_t times)
    {
      touched.clear ();
      marking.apply (net.effect (t), [this] (size_t p, T oldv, T newv) {
        enabled.onPlaceChanged (p, oldv, newv);
        touched.push_back ({ p, newv > oldv });
      }, static_cast<long long> (times));
      if (recordTrace) trace.insert (trace.end (), times, t);
      if (tracker) tracker->onFired (t, times);
    }

    /** Evaluate an open target here; record a bound's value, claim a target that holds. */
    void check (uint32_t id, WalkResult &result)
    {
      if (targets.isSolved (id)) return;
      ++result.stats.targetChecks;
      const Target<T> &tg = targets.target (id);
      if (tg.isBound ()) {
        long long v = tg.value (marking);
        if (v <= published[id]) return;
        published[id] = v;
        targets.recordValue (id, v);
        if (!tg.hasLimit () || v < tg.getLimit ()) return;
      } else if (!tg.reached (marking, enabled)) {
        return;
      }
      if (!targets.claim (id)) return;
      ++result.claims;
      if (id == focus) result.found = true;
      if (onClaim) onClaim (id, marking, recordTrace && !fromPool ? &trace : nullptr);
    }

    /** Every own open target, dropping the solved ones met on the way. */
    void checkAll (WalkResult &result)
    {
      for (size_t i = 0; i < own.size ();) {
        if (targets.isSolved (own[i])) { dropAt (own, i); continue; }
        check (own[i], result);
        ++i;
      }
    }

    /**
     * The own open targets a place changed by the last firing can have made
     * hold: those waiting for an increase of a place that grew, for a decrease
     * of one that shrank. Solved ids met on the way leave the lists.
     */
    void checkTouched (WalkResult &result)
    {
      ++epoch;
      candidates.clear ();
      for (const auto &pc : touched) {
        std::vector<uint32_t> &list = pc.second ? up[pc.first] : down[pc.first];
        for (size_t i = 0; i < list.size ();) {
          uint32_t id = list[i];
          if (targets.isSolved (id)) { dropAt (list, i); continue; }
          ++i;
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
    /**
     * A walker over all targets, or over `subset` when given (the sweep over a
     * large set is split between threads; a walker whose subset is all solved
     * takes over the open targets).
     */
    Walker (const WalkNet<T> &n, TargetSet<T> &tgs, uint32_t focusTarget, Strategy<T> &st, uint64_t seed,
            const std::vector<uint32_t> *subset = nullptr)
        : net (n), targets (tgs), focus (focusTarget), strategy (st), rng (seed),
          initialMarking (n.initialMarking ()), initialEnabled (n),
          marking (n.initialMarking ()), enabled (n),
          published (tgs.size (), std::numeric_limits<long long>::min ()), stamp (tgs.size (), 0),
          up (tgs.places ()), down (tgs.places ())
    {
      initialEnabled.initialize (initialMarking);
      enabled.assign (initialEnabled);
      if (subset) {
        buildIndex (*subset);
      } else {
        std::vector<uint32_t> all (tgs.size ());
        for (uint32_t i = 0; i < all.size (); ++i) all[i] = i;
        buildIndex (all);
      }
    }

    void setOnClaim (OnClaim cb)
    {
      onClaim = std::move (cb);
    }
    void setPool (SharedPool<T> *p)
    {
      pool = p;
    }
    void setSaturate (bool s)
    {
      saturate = s;
    }
    void setTracker (NoveltyTracker<T> *t)
    {
      tracker = t;
    }
    void setKnowledge (const Knowledge *k)
    {
      knowledge = k;
    }
    /** Ends runs; without one, the step budget of the WalkBudget applies. */
    void setRestartPolicy (const RestartPolicy *p)
    {
      restartPolicy = p;
    }

    /**
     * The walk as a resumable loop: begin() takes the budget and the initial
     * marking, runSlice() advances by at most maxSteps steps or until the
     * slice deadline and says whether the walk is over, finish() closes the
     * statistics. run() is the three in a row.
     */
    void begin (const WalkBudget &b)
    {
      budget = b;
      result = WalkResult ();
      recordTrace = budget.recordTrace;
      fallback = StepBudget (budget.runLength);
      policy = restartPolicy ? restartPolicy : &fallback;
      composite = dynamic_cast<const AnyOf*> (policy);
      runSteps = 0;
      iterations = 0;
      activeMs = 0;
      runStartMs = 0;
      over = false;
      fromPool = false;
      marking.assign (initialMarking);
      enabled.assign (initialEnabled);
      trace.clear ();
      strategy.onReset ();
      checkAll (result);
    }

    /** Advance the walk; true when it is over (targets done, budget spent, stop flag). */
    bool runSlice (uint64_t maxSteps, std::chrono::steady_clock::time_point sliceDeadline, SliceReport *report = nullptr)
    {
      using clock = std::chrono::steady_clock;
      auto sliceStart = clock::now ();
      WalkStats &st = result.stats;
      uint64_t stepsAtStart = st.steps, claimsAtStart = result.claims;
      uint64_t ms = activeMs; // running time, refreshed every 64 iterations
      WalkContext<T> ctx { net, marking, enabled, rng, knowledge, tracker ? &tracker->localCounts () : nullptr };
      bool capped = false;
      while (!over) {
        if (done ()) { over = true; break; }
        if (enabled.empty ()) {
          checkDeadlocks (result);
          if (done ()) { over = true; break; }
          ++st.deadEnds;
        }
        RunView view { runSteps, ms - runStartMs, tracker ? tracker->stepsSinceNovelty () : 0 };
        bool policyEnd = !enabled.empty () && policy->shouldRestart (view);
        if (enabled.empty () || policyEnd) {
          ++st.resets;
          if (policyEnd) ++st.policyResets;
          runSteps = 0;
          runStartMs = ms;
          bool inPlace = policyEnd && (composite ? composite->keepsMarking (view) : policy->keepsMarking ());
          if (inPlace) {
            ++st.inPlaceResets;
            resetInPlace ();
            refill ();
            continue;
          }
          reset (&st);
          refill ();
          // back at the initial marking nothing has changed since the first check
          if (fromPool) checkAll (result);
          continue;
        }
        if (st.steps - stepsAtStart >= maxSteps) break;
        // the clock every 64 steps: a step costs microseconds on most nets, milliseconds on hub-dense ones
        if ((iterations++ & 63) == 0) {
          if (budget.stop && budget.stop->load (std::memory_order_relaxed)) { over = true; break; }
          if (budget.maxSteps && st.steps >= budget.maxSteps) { over = true; break; }
          auto now = clock::now ();
          ms = activeMs + static_cast<uint64_t> (std::chrono::duration_cast<std::chrono::milliseconds> (now - sliceStart).count ());
          if (budget.timeoutMillis && ms >= budget.timeoutMillis) { over = true; break; }
          if (now >= sliceDeadline) { capped = true; break; }
          if (tracker && (iterations & 4095) == 1) tracker->merge ();
        }
        uint32_t t = strategy.choose (ctx);
        if (ctx.badStart) {
          // the strategy found the run's start hopeless: a pooled state leaves the pool, and the walk restarts
          ctx.badStart = false;
          if (fromPool && pool) {
            pool->evict (poolEntry.id);
            ++st.condemned;
            t = RESTART;
          }
        }
        if (t == RESTART) {
          ++st.stalls;
          ++st.resets;
          runSteps = 0;
          runStartMs = ms;
          reset (&st);
          refill ();
          if (fromPool) checkAll (result);
          continue;
        }
        uint64_t times = saturate ? maxFirings (t) : 1;
        fire (t, times);
        if (times > 1) ++st.saturations;
        strategy.onFired (t, times);
        st.steps += times;
        runSteps += times;
        checkTouched (result);
        // a dead marking is no restart point: nothing can be explored from it
        if (tracker && !enabled.empty ()) tracker->observe (marking);
      }
      uint64_t micros = static_cast<uint64_t> (std::chrono::duration_cast<std::chrono::microseconds> (clock::now () - sliceStart).count ());
      activeMs += micros / 1000;
      if (report) {
        report->steps = st.steps - stepsAtStart;
        report->claims = result.claims - claimsAtStart;
        report->micros = micros;
        report->capped = capped && !over;
        report->finished = over;
      }
      return over;
    }

    WalkResult finish ()
    {
      WalkStats &st = result.stats;
      publishRun ();
      if (tracker) {
        tracker->merge ();
        st.distinctFired = tracker->distinctFired ();
        st.rareEvents = tracker->rareEventCount ();
      }
      st.millis = activeMs;
      st.arcVisits = enabled.arcVisits;
      st.flips = enabled.flips;
      return result;
    }

    WalkResult run (const WalkBudget &b)
    {
      begin (b);
      runSlice (std::numeric_limits<uint64_t>::max (), std::chrono::steady_clock::time_point::max ());
      return finish ();
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
        fire (t, 1);
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
