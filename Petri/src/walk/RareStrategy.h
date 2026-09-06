/*
 * RareStrategy.h
 *
 * Transition choice by rarity and age: among a few enabled transitions drawn
 * uniformly, fire the one fired the fewest times by any thread, the longest
 * enabled among equals. See algorithm.md.
 */
#ifndef PETRI_WALK_RARESTRATEGY_H_
#define PETRI_WALK_RARESTRATEGY_H_

#include <algorithm>
#include <cstdint>

#include "walk/Strategy.h"

namespace petri::walk
{

/**
 * Rarity and age: draw a few enabled transitions uniformly and fire the one
 * fired the fewest times by anyone (the walker's own count added to the
 * shared one), the longest enabled among equals. Uniform choice makes a
 * transition that needs a rare conjunction exponentially rare; this makes it
 * the preferred move the moment it is enabled, at O(sample) per step and
 * without scanning an enabled set that may hold tens of thousands.
 */
template<typename T>
  class RareStrategy : public Strategy<T>
  {
    size_t sample;

  public:
    explicit RareStrategy (size_t k)
        : sample (k == 0 ? 8 : k)
    {
    }

    uint32_t choose (WalkContext<T> &ctx) override
    {
      size_t n = ctx.enabled.size ();
      size_t k = std::min (sample, n);
      uint32_t best = ctx.enabled.at (ctx.rng () % n);
      uint64_t bestFired = fired (ctx, best);
      uint32_t bestAge = ctx.enabled.age (best);
      for (size_t i = 1; i < k; ++i) {
        uint32_t t = ctx.enabled.at (ctx.rng () % n);
        uint64_t f = fired (ctx, t);
        uint32_t a = ctx.enabled.age (t);
        if (f < bestFired || (f == bestFired && a > bestAge)) {
          best = t;
          bestFired = f;
          bestAge = a;
        }
      }
      return best;
    }

  private:
    static uint64_t fired (const WalkContext<T> &ctx, uint32_t t)
    {
      uint64_t f = ctx.knowledge ? ctx.knowledge->fired (t) : 0;
      if (ctx.firedLocal) f += (*ctx.firedLocal)[t];
      return f;
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_RARESTRATEGY_H_ */
