/*
 * Strategy.h
 *
 * Choice of the next transition to fire among the enabled ones.
 */
#ifndef PETRI_WALK_STRATEGY_H_
#define PETRI_WALK_STRATEGY_H_

#include <cstdint>
#include <limits>
#include <random>

#include "walk/EnabledSet.h"
#include "walk/Marking.h"
#include "walk/WalkNet.h"

namespace petri::walk
{

template<typename T>
  struct WalkContext
  {
    const WalkNet<T> &net;
    const Marking<T> &marking;
    const EnabledSet<T> &enabled;
    std::mt19937_64 &rng;
  };

/** Returned by Strategy::choose to ask the walker to restart the run. */
constexpr uint32_t RESTART = std::numeric_limits<uint32_t>::max ();

template<typename T>
  class Strategy
  {
  public:
    virtual ~Strategy () = default;
    /** Called with a non-empty enabled set; returns a transition of it, or RESTART. */
    virtual uint32_t choose (WalkContext<T> &ctx) = 0;
    virtual void onReset ()
    {
    }
    /** The walker fired t, times times in one update (more than once under saturation). */
    virtual void onFired (uint32_t, uint64_t)
    {
    }
    /**
     * Best state of the current run according to the strategy's heuristic,
     * if it has one: fills m and h and returns true. Random strategies
     * return false.
     */
    virtual bool bestOfRun (SparseArray<T> &, uint64_t &) const
    {
      return false;
    }
  };

/** Uniform choice among enabled transitions. */
template<typename T>
  class RandomStrategy : public Strategy<T>
  {
  public:
    uint32_t choose (WalkContext<T> &ctx) override
    {
      size_t n = ctx.enabled.size ();
      return ctx.enabled.at (ctx.rng () % n);
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_STRATEGY_H_ */
