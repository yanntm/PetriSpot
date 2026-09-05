/*
 * DeadlockStrategy.h
 *
 * Greedy choice on the number of enabled transitions of the successor
 * marking: the distance to a dead marking is the size of the enabled set, and
 * zero is exactly the goal. Random tie-breaking, an epsilon share of uniform
 * moves and a stall restart, as in BestFirstStrategy.
 *
 * The successor count is obtained without building the successor marking. A
 * transition changes status only when one of the places the candidate touches
 * crosses one of its consumer weights, so the count moves by the flips over
 * the consumers of those places. A transition can be touched through several
 * places at once, so the per transition deltas are accumulated before the
 * flips are read, against the unsatisfied-arc counters the enabled set keeps.
 */
#ifndef PETRI_WALK_DEADLOCKSTRATEGY_H_
#define PETRI_WALK_DEADLOCKSTRATEGY_H_

#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

#include "walk/Strategy.h"

namespace petri::walk
{

template<typename T>
  class DeadlockStrategy : public Strategy<T>
  {
    static constexpr uint64_t NO_COUNT = std::numeric_limits<uint64_t>::max ();

    unsigned epsilonPercent; // share of uniformly random moves, in percent
    size_t sampleSize;       // at most this many candidates scored per step (0: all)
    uint64_t stallLimit;     // steps without improvement before a restart (0: never)
    uint64_t sinceImprovement = 0;
    std::vector<uint32_t> best;

    // Per candidate scratch. Epoch stamps keep it clear of the previous
    // candidate without touching anything but the transitions it reached.
    std::vector<int32_t> delta;
    std::vector<uint32_t> stamp;
    std::vector<uint32_t> touched;
    uint32_t epoch = 0;

  public:
    /** Instrumentation: fewest enabled transitions seen, overall and per run. */
    uint64_t minEnabled = NO_COUNT;
    uint64_t minEnabledThisRun = NO_COUNT;
    SparseArray<T> bestMarkingThisRun;

    DeadlockStrategy (const WalkNet<T> &net, unsigned epsilon = 10,
                      size_t sample = 0, uint64_t stall = 0)
        : epsilonPercent (epsilon), sampleSize (sample), stallLimit (stall),
          delta (net.transitionCount (), 0), stamp (net.transitionCount (), 0)
    {
    }

    void onReset () override
    {
      minEnabledThisRun = NO_COUNT;
      sinceImprovement = 0;
    }

    bool bestOfRun (SparseArray<T> &m, uint64_t &h) const override
    {
      if (minEnabledThisRun == NO_COUNT) return false;
      m = bestMarkingThisRun;
      h = minEnabledThisRun;
      return true;
    }

    /** Enabled transitions of the marking reached by firing t. */
    uint64_t successorCount (WalkContext<T> &ctx, uint32_t t)
    {
      ++epoch;
      touched.clear ();
      const SparseArray<T> &eff = ctx.net.effect (t);
      for (size_t i = 0, e = eff.size (); i < e; ++i) {
        const T d = eff.valueAt (i);
        if (d == 0) continue;
        const size_t p = eff.keyAt (i);
        const auto &cs = ctx.net.consumersOf (p);
        if (cs.empty ()) continue;
        const T oldv = ctx.marking.get (p);
        const T newv = static_cast<T> (oldv + d);
        const T maxW = ctx.net.maxWeight (p);
        if (oldv >= maxW && newv >= maxW) continue;
        T lo, hi;
        bool becomeSat;
        if (newv > oldv) {
          lo = oldv;
          hi = newv;
          becomeSat = true;
        } else {
          lo = newv;
          hi = oldv;
          becomeSat = false;
        }
        auto first = std::upper_bound (cs.begin (), cs.end (), lo,
            [] (T v, const typename WalkNet<T>::Consumer &c) { return v < c.weight; });
        auto last = std::upper_bound (first, cs.end (), hi,
            [] (T v, const typename WalkNet<T>::Consumer &c) { return v < c.weight; });
        for (auto it = first; it != last; ++it) {
          const uint32_t c = it->transition;
          if (stamp[c] != epoch) {
            stamp[c] = epoch;
            delta[c] = 0;
            touched.push_back (c);
          }
          delta[c] += becomeSat ? -1 : 1;
        }
      }
      long long count = static_cast<long long> (ctx.enabled.size ());
      for (uint32_t c : touched) {
        const long long before = static_cast<long long> (ctx.enabled.unsatOf (c));
        const long long after = before + delta[c];
        const bool was = (before == 0);
        const bool now = (after <= 0);
        if (was && !now) --count;
        else if (!was && now) ++count;
      }
      return count < 0 ? 0 : static_cast<uint64_t> (count);
    }

    uint32_t choose (WalkContext<T> &ctx) override
    {
      const size_t n = ctx.enabled.size ();
      if (stallLimit && sinceImprovement >= stallLimit) return RESTART;
      if (n == 1 || (epsilonPercent > 0 && ctx.rng () % 100 < epsilonPercent)) {
        ++sinceImprovement;
        return ctx.enabled.at (ctx.rng () % n);
      }
      best.clear ();
      uint64_t bestCount = NO_COUNT;
      const size_t count = sampleSize == 0 || sampleSize >= n ? n : sampleSize;
      for (size_t i = 0; i < count; ++i) {
        const uint32_t t = count == n ? ctx.enabled.at (i) : ctx.enabled.at (ctx.rng () % n);
        const uint64_t d = successorCount (ctx, t);
        if (d < bestCount) {
          bestCount = d;
          best.clear ();
          best.push_back (t);
        } else if (d == bestCount) {
          best.push_back (t);
        }
      }
      ++sinceImprovement;
      if (bestCount < minEnabledThisRun) {
        minEnabledThisRun = bestCount;
        sinceImprovement = 0;
        bestMarkingThisRun = ctx.marking.sparse ();
        if (bestCount < minEnabled) minEnabled = bestCount;
      }
      return best[ctx.rng () % best.size ()];
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_DEADLOCKSTRATEGY_H_ */
