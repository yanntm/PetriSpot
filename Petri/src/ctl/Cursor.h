/*
 * Cursor.h
 *
 * A marking with its enabled set: the state a solve works on. Copied when a
 * sub-obligation is opened (the copy costs the transition count), fired and
 * reverted in place along a walk or a DFS.
 */
#ifndef PETRI_CTL_CURSOR_H_
#define PETRI_CTL_CURSOR_H_

#include <cstdint>

#include "walk/EnabledSet.h"
#include "walk/Marking.h"
#include "walk/WalkNet.h"

namespace petri::ctl
{

template<typename T>
  struct Cursor
  {
    const petri::walk::WalkNet<T> *net;
    petri::walk::Marking<T> marking;
    petri::walk::EnabledSet<T> enabled;

    explicit Cursor (const petri::walk::WalkNet<T> &n)
        : net (&n), marking (n.initialMarking ()), enabled (n)
    {
      enabled.initialize (marking);
    }
    Cursor (const Cursor &) = default;
    Cursor& operator= (const Cursor &) = default;

    bool deadlock () const
    {
      return enabled.empty ();
    }

    void fire (uint32_t t)
    {
      marking.apply (net->effect (t), [this] (size_t p, T o, T n) { enabled.onPlaceChanged (p, o, n); }, 1);
    }

    /** Undo a fire of t (the marking is the one fire produced). */
    void revert (uint32_t t)
    {
      marking.apply (net->effect (t), [this] (size_t p, T o, T n) { enabled.onPlaceChanged (p, o, n); }, -1);
    }

    const SparseArray<T>& state () const
    {
      return marking.sparse ();
    }
  };

} // namespace petri::ctl

#endif /* PETRI_CTL_CURSOR_H_ */
