/*
 * ParikhStrategy.h
 *
 * Hint-driven choice: fire transitions whose effect class still has firings
 * left in a Parikh vector; restart when none is enabled; relax the
 * restriction with a probability growing with the number of restarts. See
 * algorithm.md.
 */
#ifndef PETRI_WALK_PARIKHSTRATEGY_H_
#define PETRI_WALK_PARIKHSTRATEGY_H_

#include <algorithm>
#include <cstdint>
#include <vector>

#include "expr/Hint.h"
#include "walk/Strategy.h"
#include "walk/WalkNet.h"

namespace petri::walk
{

template<typename T>
  class ParikhStrategy : public Strategy<T>
  {
    const WalkNet<T> &net;
    std::vector<long long> initial;   // per effect class: firings suggested
    std::vector<long long> remaining; // per effect class: firings left in this run
    uint64_t resets = 0;
    std::vector<uint32_t> candidates;

  public:
    /**
     * Counts are given on transitions; a count lands on the transition's
     * effect class (the largest when several members are named), since the
     * state equation does not distinguish transitions with the same effect.
     */
    ParikhStrategy (const WalkNet<T> &n, const petri::expr::ParikhHint &hint)
        : net (n), initial (n.effectClassCount (), 0), remaining (n.effectClassCount (), 0)
    {
      for (const auto &c : hint.counts) {
        long long &slot = initial[net.effectClassOf (c.first)];
        slot = std::max (slot, c.second);
      }
      remaining = initial;
    }

    void onReset () override
    {
      remaining = initial;
      ++resets;
    }

    void onFired (uint32_t t, uint64_t times) override
    {
      long long &left = remaining[net.effectClassOf (t)];
      left = std::max (0LL, left - static_cast<long long> (times));
    }

    uint32_t choose (WalkContext<T> &ctx) override
    {
      const size_t n = ctx.enabled.size ();
      // the Java decay: skip the restriction with probability resets per mille
      bool relaxed = (ctx.rng () % 1000) < std::min<uint64_t> (resets, 1000);
      if (relaxed) return ctx.enabled.at (ctx.rng () % n);
      candidates.clear ();
      for (size_t i = 0; i < n; ++i) {
        uint32_t t = ctx.enabled.at (i);
        if (remaining[net.effectClassOf (t)] > 0) candidates.push_back (t);
      }
      if (candidates.empty ()) return RESTART;
      return candidates[ctx.rng () % candidates.size ()];
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_PARIKHSTRATEGY_H_ */
