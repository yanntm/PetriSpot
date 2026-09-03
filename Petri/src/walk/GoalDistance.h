/*
 * GoalDistance.h
 *
 * Estimated distance of a marking to the goal of a walk; zero iff the goal
 * holds. MarkingDistance evaluates the TAPAAL-style expression distance.
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

} // namespace petri::walk

#endif /* PETRI_WALK_GOALDISTANCE_H_ */
