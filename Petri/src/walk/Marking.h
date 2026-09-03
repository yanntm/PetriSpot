/*
 * Marking.h
 *
 * Sparse marking of a net, mutated in place by transition effects.
 * See algorithm.md.
 */
#ifndef PETRI_WALK_MARKING_H_
#define PETRI_WALK_MARKING_H_

#include <vector>

#include "core/SparseArray.h"

namespace petri::walk
{

template<typename T>
  class Marking
  {
    SparseArray<T> m;

  public:
    Marking () = default;
    explicit Marking (const std::vector<T> &dense)
        : m (dense)
    {
    }

    T get (size_t p) const
    {
      return m.get (p);
    }
    const SparseArray<T>& sparse () const
    {
      return m;
    }
    /** Copy; reuses this marking's storage when large enough. */
    void assign (const Marking &o)
    {
      m = o.m;
    }
    bool covers (const SparseArray<T> &pre) const
    {
      return SparseArray<T>::greaterOrEqual (m, pre);
    }

    /**
     * Add effect to this marking; onChange(place, oldValue, newValue) is
     * called for every touched place, after the marking is updated.
     */
    template<class OnChange>
      void apply (const SparseArray<T> &effect, OnChange &&onChange)
      {
        for (size_t i = 0, e = effect.size (); i < e; ++i) {
          size_t p = effect.keyAt (i);
          T d = effect.valueAt (i);
          ssize_t idx = m.indexOfKey (p);
          T oldv = 0;
          if (idx >= 0) {
            oldv = m.valueAt (static_cast<size_t> (idx));
            T newv = oldv + d;
            if (newv == 0) m.removeAt (static_cast<size_t> (idx));
            else m.setValueAt (static_cast<size_t> (idx), newv);
            onChange (p, oldv, newv);
          } else {
            m.put (p, d);
            onChange (p, oldv, d);
          }
        }
      }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_MARKING_H_ */
