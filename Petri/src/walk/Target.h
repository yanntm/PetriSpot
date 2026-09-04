/*
 * Target.h
 *
 * The goal of a walk: a state predicate in normal form, a deadlock, or a
 * linear form to maximise (a bound), with an optional limit that ends it.
 */
#ifndef PETRI_WALK_TARGET_H_
#define PETRI_WALK_TARGET_H_

#include <limits>

#include "expr/Expression.h"
#include "walk/EnabledSet.h"
#include "walk/Marking.h"

namespace petri::walk
{

template<typename T>
  class Target
  {
  public:
    enum class Kind
    {
      Goal, Deadlock, Bound
    };

  private:
    Kind kind;
    petri::expr::Expression goal;    // Goal: the predicate; Bound with a limit: form >= limit
    petri::expr::LinearAtom form;    // Bound: the terms to maximise
    long long limit = std::numeric_limits<long long>::max (); // Bound: value that ends the target

  public:
    /** Reach a state satisfying g. */
    explicit Target (petri::expr::Expression g)
        : kind (Kind::Goal), goal (std::move (g))
    {
    }
    /** Reach a deadlock. */
    static Target deadlockTarget ()
    {
      Target t (petri::expr::Expression::constant (false));
      t.kind = Kind::Deadlock;
      return t;
    }
    /**
     * Maximise form; a non-negative limit is a known upper bound whose value,
     * once reached, ends the target (goal is then "form >= limit").
     */
    static Target boundTarget (petri::expr::LinearAtom f, long long knownLimit)
    {
      Target t (petri::expr::Expression::constant (false));
      t.kind = Kind::Bound;
      t.form = std::move (f);
      if (knownLimit >= 0) {
        t.limit = knownLimit;
        petri::expr::LinearAtom a = t.form;
        a.op = petri::expr::Cmp::GE;
        a.constant = knownLimit;
        t.goal = petri::expr::Expression::makeAtom (std::move (a));
      }
      return t;
    }

    Kind getKind () const
    {
      return kind;
    }
    bool isDeadlock () const
    {
      return kind == Kind::Deadlock;
    }
    bool isBound () const
    {
      return kind == Kind::Bound;
    }
    /** Bound: whether a limit is known. */
    bool hasLimit () const
    {
      return limit != std::numeric_limits<long long>::max ();
    }
    long long getLimit () const
    {
      return limit;
    }
    /** Goal, or Bound with a limit: the predicate; False for the other kinds. */
    const petri::expr::Expression& expression () const
    {
      return goal;
    }
    const petri::expr::LinearAtom& boundForm () const
    {
      return form;
    }
    /** Bound: the value of the form in m. */
    long long value (const Marking<T> &m) const
    {
      return form.value (m);
    }

    bool reached (const Marking<T> &m, const EnabledSet<T> &enabled) const
    {
      switch (kind) {
      case Kind::Deadlock: return enabled.empty ();
      case Kind::Bound: return hasLimit () && value (m) >= limit;
      case Kind::Goal: return goal.eval (m);
      }
      return false;
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_TARGET_H_ */
