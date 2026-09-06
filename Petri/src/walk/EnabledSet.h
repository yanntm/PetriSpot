/*
 * EnabledSet.h
 *
 * Enabled transitions maintained by delta: per-transition count of
 * unsatisfied input arcs and a packed list of the enabled transitions.
 * See algorithm.md.
 */
#ifndef PETRI_WALK_ENABLEDSET_H_
#define PETRI_WALK_ENABLEDSET_H_

#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

#include "walk/Marking.h"
#include "walk/WalkNet.h"

namespace petri::walk
{

template<typename T>
  class EnabledSet
  {
    static constexpr uint32_t NONE = std::numeric_limits<uint32_t>::max ();

    const WalkNet<T> *net;
    std::vector<uint32_t> unsat;
    std::vector<uint32_t> position;
    std::vector<uint32_t> list;
    std::vector<uint32_t> since;    // per transition: the clock at which it last became enabled
    uint32_t clock = 0;             // advanced at every enabling

  public:
    /** Instrumentation: consumer arcs visited and enabled-status flips. */
    uint64_t arcVisits = 0;
    uint64_t flips = 0;

  private:
    void add (uint32_t t)
    {
      position[t] = static_cast<uint32_t> (list.size ());
      list.push_back (t);
      since[t] = ++clock;
    }
    void remove (uint32_t t)
    {
      uint32_t pos = position[t];
      uint32_t last = list.back ();
      list[pos] = last;
      position[last] = pos;
      list.pop_back ();
      position[t] = NONE;
    }

  public:
    explicit EnabledSet (const WalkNet<T> &n)
        : net (&n), unsat (n.transitionCount (), 0),
          position (n.transitionCount (), NONE), since (n.transitionCount (), 0)
    {
      list.reserve (n.transitionCount ());
    }

    /** Recompute everything from the marking: O(arcs). */
    void initialize (const Marking<T> &m)
    {
      list.clear ();
      std::fill (position.begin (), position.end (), NONE);
      clock = 0;
      const size_t nbT = net->transitionCount ();
      for (size_t t = 0; t < nbT; ++t) {
        const SparseArray<T> &pre = net->pre (t);
        uint32_t missing = 0;
        for (size_t i = 0; i < pre.size (); ++i) {
          if (m.get (pre.keyAt (i)) < pre.valueAt (i)) ++missing;
        }
        unsat[t] = missing;
        if (missing == 0) add (static_cast<uint32_t> (t));
      }
    }

    /** Copy; reuses storage. */
    void assign (const EnabledSet &o)
    {
      net = o.net;
      unsat = o.unsat;
      position = o.position;
      list = o.list;
      since = o.since;
      clock = o.clock;
    }

    /** Delta update after place p moved from oldv to newv. */
    void onPlaceChanged (size_t p, T oldv, T newv)
    {
      const auto &cs = net->consumersOf (p);
      if (cs.empty ()) return;
      T maxW = net->maxWeight (p);
      if (oldv >= maxW && newv >= maxW) return;
      T lo, hi; // consumers with weight in (lo, hi] flip
      bool becomeSat;
      if (newv > oldv) {
        lo = oldv;
        hi = newv;
        becomeSat = true;
      } else {
        lo = newv;
        hi = oldv;
        becomeSat = false;
      }
      auto first = std::upper_bound (cs.begin (), cs.end (), lo,
          [] (T v, const typename WalkNet<T>::Consumer &c) { return v < c.weight; });
      auto last = std::upper_bound (first, cs.end (), hi,
          [] (T v, const typename WalkNet<T>::Consumer &c) { return v < c.weight; });
      arcVisits += static_cast<uint64_t> (last - first);
      for (auto it = first; it != last; ++it) {
        uint32_t t = it->transition;
        if (becomeSat) {
          if (--unsat[t] == 0) { add (t); ++flips; }
        } else {
          if (unsat[t]++ == 0) { remove (t); ++flips; }
        }
      }
    }

    /** Unsatisfied pre-arcs of t; t is enabled exactly when this is zero. */
    uint32_t unsatOf (size_t t) const
    {
      return unsat[t];
    }

    size_t size () const
    {
      return list.size ();
    }
    bool empty () const
    {
      return list.empty ();
    }
    uint32_t at (size_t i) const
    {
      return list[i];
    }
    bool isEnabled (size_t t) const
    {
      return unsat[t] == 0;
    }
    /** Enablings that happened since t last became enabled: the larger, the longer t has waited. */
    uint32_t age (size_t t) const
    {
      return clock - since[t];
    }
    const std::vector<uint32_t>& transitions () const
    {
      return list;
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_ENABLEDSET_H_ */
