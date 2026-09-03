/*
 * BestFirstStrategy.h
 *
 * Greedy choice on the estimated GoalDistance of the successor marking, with random tie-breaking and an epsilon share of uniform
 * random moves.
 */
#ifndef PETRI_WALK_BESTFIRSTSTRATEGY_H_
#define PETRI_WALK_BESTFIRSTSTRATEGY_H_

#include <cstdint>
#include <vector>

#include "expr/Distance.h"
#include "walk/GoalDistance.h"
#include "walk/Strategy.h"

namespace petri::walk
{

template<typename T>
  class BestFirstStrategy : public Strategy<T>
  {
    const GoalDistance<T> &goal;
    unsigned epsilonPercent; // share of uniformly random moves, in percent
    size_t sampleSize;       // at most this many candidates scored per step (0: all)
    uint64_t stallLimit;     // steps without improvement before a restart (0: never)
    uint64_t sinceImprovement = 0;
    std::vector<uint32_t> best;

  public:
    /** Instrumentation: smallest successor distance chosen, overall and per run. */
    uint64_t minDistance = petri::expr::INFINITE_DISTANCE;
    uint64_t minDistanceThisRun = petri::expr::INFINITE_DISTANCE;
    uint64_t runsReachingMin = 0;
    SparseArray<T> bestMarking; // marking from which the best successor was chosen

    void onReset () override
    {
      minDistanceThisRun = petri::expr::INFINITE_DISTANCE;
      sinceImprovement = 0;
    }

    BestFirstStrategy (const GoalDistance<T> &g, unsigned epsilon = 10,
                       size_t sample = 0, uint64_t stall = 0)
        : goal (g), epsilonPercent (epsilon), sampleSize (sample), stallLimit (stall)
    {
    }

    uint32_t choose (WalkContext<T> &ctx) override
    {
      const size_t n = ctx.enabled.size ();
      if (stallLimit && sinceImprovement >= stallLimit) return RESTART;
      if (n == 1 || (epsilonPercent > 0 && ctx.rng () % 100 < epsilonPercent)) {
        ++sinceImprovement;
        return ctx.enabled.at (ctx.rng () % n);
      }
      Marking<T> &m = const_cast<Marking<T>&> (ctx.marking); // restored by peek
      best.clear ();
      uint64_t bestDist = petri::expr::INFINITE_DISTANCE + 1;
      size_t count = sampleSize == 0 || sampleSize >= n ? n : sampleSize;
      for (size_t i = 0; i < count; ++i) {
        uint32_t t = count == n ? ctx.enabled.at (i) : ctx.enabled.at (ctx.rng () % n);
        uint64_t d = m.peek (ctx.net.effect (t), [this] (const Marking<T> &mm) {
          return goal.of (mm);
        });
        if (d < bestDist) {
          bestDist = d;
          best.clear ();
          best.push_back (t);
        } else if (d == bestDist) {
          best.push_back (t);
        }
      }
      ++sinceImprovement;
      if (bestDist < minDistanceThisRun) {
        minDistanceThisRun = bestDist;
        sinceImprovement = 0;
        if (bestDist < minDistance) {
          minDistance = bestDist;
          runsReachingMin = 1;
          bestMarking = m.sparse ();
        } else if (bestDist == minDistance) {
          ++runsReachingMin;
        }
      }
      return best[ctx.rng () % best.size ()];
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_BESTFIRSTSTRATEGY_H_ */
