/*
 * Target.h
 *
 * The goal of a walk: a state predicate in normal form, or deadlock.
 */
#ifndef PETRI_WALK_TARGET_H_
#define PETRI_WALK_TARGET_H_

#include "expr/Expression.h"
#include "walk/EnabledSet.h"
#include "walk/Marking.h"

namespace petri::walk
{

template<typename T>
  class Target
  {
    petri::expr::Expression goal;
    bool deadlock;

  public:
    /** Reach a state satisfying goal. */
    explicit Target (petri::expr::Expression g)
        : goal (std::move (g)), deadlock (false)
    {
    }
    /** Reach a deadlock. */
    static Target deadlockTarget ()
    {
      Target t (petri::expr::Expression::constant (false));
      t.deadlock = true;
      return t;
    }

    bool isDeadlock () const
    {
      return deadlock;
    }
    const petri::expr::Expression& expression () const
    {
      return goal;
    }

    bool reached (const Marking<T> &m, const EnabledSet<T> &enabled) const
    {
      if (deadlock) return enabled.empty ();
      return goal.eval (m);
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_TARGET_H_ */
