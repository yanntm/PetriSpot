/*
 * GoalDistance.h
 *
 * Estimated distance of a marking to the goal of a walk; zero iff the goal
 * holds. MarkingDistance evaluates the TAPAAL-style expression distance;
 * BoundDistance ranks markings by the value of a linear form to maximise
 * (never zero: a bound without a known limit has no goal to reach).
 */
#ifndef PETRI_WALK_GOALDISTANCE_H_
#define PETRI_WALK_GOALDISTANCE_H_

#include <cstdint>

#include "expr/Distance.h"
#include "expr/Expression.h"
#include "walk/Marking.h"

namespace petri::walk
{

template<typename T>
  class GoalDistance
  {
  public:
    virtual ~GoalDistance () = default;
    virtual uint64_t of (const Marking<T> &m) const = 0;
  };

template<typename T>
  class MarkingDistance : public GoalDistance<T>
  {
    const petri::expr::Expression &goal;

  public:
    explicit MarkingDistance (const petri::expr::Expression &g)
        : goal (g)
    {
    }
    uint64_t of (const Marking<T> &m) const override
    {
      return petri::expr::distance (goal, m);
    }
  };

/**
 * For a bound target without a limit: OFFSET minus the value of the form, so
 * that a larger value is a smaller distance and improvements are monotone
 * (which the stall detection of best-first and the shared pool rely on).
 */
template<typename T>
  class BoundDistance : public GoalDistance<T>
  {
    const petri::expr::LinearAtom &form;
    static constexpr uint64_t OFFSET = petri::expr::INFINITE_DISTANCE / 2;

  public:
    explicit BoundDistance (const petri::expr::LinearAtom &f)
        : form (f)
    {
    }
    uint64_t of (const Marking<T> &m) const override
    {
      long long v = form.value (m);
      if (v <= 0) return OFFSET + static_cast<uint64_t> (-v);
      return OFFSET - static_cast<uint64_t> (v);
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_GOALDISTANCE_H_ */
