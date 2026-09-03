/*
 * BestFirstStrategy.h
 *
 * Greedy choice on the estimated distance of the successor marking to the
 * goal expression, with random tie-breaking and an epsilon share of uniform
 * random moves.
 */
#ifndef PETRI_WALK_BESTFIRSTSTRATEGY_H_
#define PETRI_WALK_BESTFIRSTSTRATEGY_H_

#include <cstdint>
#include <vector>

#include "expr/Distance.h"
#include "expr/Expression.h"
#include "walk/Strategy.h"

namespace petri::walk
{

template<typename T>
  class BestFirstStrategy : public Strategy<T>
  {
    const petri::expr::Expression &goal;
    unsigned epsilonPercent; // share of uniformly random moves, in percent
    size_t sampleSize;       // at most this many candidates scored per step (0: all)
    std::vector<uint32_t> best;

  public:
    /** Instrumentation: smallest successor distance chosen, overall and per run. */
    uint64_t minDistance = petri::expr::INFINITE_DISTANCE;
    uint64_t minDistanceThisRun = petri::expr::INFINITE_DISTANCE;
    uint64_t runsReachingMin = 0;

    void onReset () override
    {
      minDistanceThisRun = petri::expr::INFINITE_DISTANCE;
    }

    BestFirstStrategy (const petri::expr::Expression &g, unsigned epsilon = 10,
                       size_t sample = 0)
        : goal (g), epsilonPercent (epsilon), sampleSize (sample)
    {
    }

    uint32_t choose (WalkContext<T> &ctx) override
    {
      const size_t n = ctx.enabled.size ();
      if (n == 1 || (epsilonPercent > 0 && ctx.rng () % 100 < epsilonPercent)) {
        return ctx.enabled.at (ctx.rng () % n);
      }
      Marking<T> &m = const_cast<Marking<T>&> (ctx.marking); // restored by peek
      best.clear ();
      uint64_t bestDist = petri::expr::INFINITE_DISTANCE + 1;
      size_t count = sampleSize == 0 || sampleSize >= n ? n : sampleSize;
      for (size_t i = 0; i < count; ++i) {
        uint32_t t = count == n ? ctx.enabled.at (i) : ctx.enabled.at (ctx.rng () % n);
        uint64_t d = m.peek (ctx.net.effect (t), [this] (const Marking<T> &mm) {
          return petri::expr::distance (goal, mm);
        });
        if (d < bestDist) {
          bestDist = d;
          best.clear ();
          best.push_back (t);
        } else if (d == bestDist) {
          best.push_back (t);
        }
      }
      if (bestDist < minDistanceThisRun) {
        minDistanceThisRun = bestDist;
        if (bestDist < minDistance) {
          minDistance = bestDist;
          runsReachingMin = 1;
        } else if (bestDist == minDistance) {
          ++runsReachingMin;
        }
      }
      return best[ctx.rng () % best.size ()];
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_BESTFIRSTSTRATEGY_H_ */
