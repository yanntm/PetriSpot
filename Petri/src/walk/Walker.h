/*
 * Walker.h
 *
 * One explicit walk toward a Target: restart loop, step budget, wall clock,
 * witness trace. See algorithm.md.
 */
#ifndef PETRI_WALK_WALKER_H_
#define PETRI_WALK_WALKER_H_

#include <chrono>
#include <cstdint>
#include <random>
#include <vector>

#include "walk/EnabledSet.h"
#include "walk/Marking.h"
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
  uint64_t millis = 0;
  uint64_t arcVisits = 0; // consumer arcs visited by the enabled-set update
  uint64_t flips = 0;     // enabled-status changes
};

struct WalkResult
{
  bool found = false;
  std::vector<uint32_t> trace; // transitions fired from the initial marking
  WalkStats stats;
};

struct WalkBudget
{
  uint64_t maxSteps = 0;       // 0: unlimited
  uint64_t runLength = 1000000; // steps between resets
  uint64_t timeoutMillis = 0;  // 0: unlimited
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

    void reset ()
    {
      marking.assign (initialMarking);
      enabled.assign (initialEnabled);
      trace.clear ();
      strategy.onReset ();
    }

    void fire (uint32_t t)
    {
      marking.apply (net.effect (t), [this] (size_t p, T oldv, T newv) {
        enabled.onPlaceChanged (p, oldv, newv);
      });
      trace.push_back (t);
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
      WalkStats &st = result.stats;
      WalkContext<T> ctx { net, marking, enabled, rng };
      uint64_t runSteps = 0;
      reset ();
      for (;;) {
        if (target.reached (marking, enabled)) {
          result.found = true;
          result.trace = trace;
          break;
        }
        if (enabled.empty () || runSteps >= budget.runLength) {
          if (enabled.empty ()) ++st.deadEnds;
          ++st.resets;
          runSteps = 0;
          reset ();
          continue;
        }
        if ((st.steps & 1023) == 0) {
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
          reset ();
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
      reset ();
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
  };

} // namespace petri::walk

#endif /* PETRI_WALK_WALKER_H_ */
