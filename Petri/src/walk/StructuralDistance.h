/*
 * StructuralDistance.h
 *
 * Goal distance combining the expression distance with a structural
 * estimate for atoms of the form "place >= k": when the place lacks tokens,
 * the estimate is the smallest number of transitions along which a token
 * present in the marking can travel to that place, computed once by a
 * backward breadth-first search from the goal place over producing
 * transitions (each pre-place of a producer is one hop further). Optimistic
 * by construction: the other pre-places of a transition are not accounted.
 */
#ifndef PETRI_WALK_STRUCTURALDISTANCE_H_
#define PETRI_WALK_STRUCTURALDISTANCE_H_

#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>

#include "core/MatrixCol.h"
#include "expr/Distance.h"
#include "expr/Expression.h"
#include "walk/GoalDistance.h"
#include "walk/Marking.h"
#include "walk/WalkNet.h"

namespace petri::walk
{

template<typename T>
  class StructuralDistance : public GoalDistance<T>
  {
    using Expression = petri::expr::Expression;
    static constexpr uint16_t UNREACHED = std::numeric_limits<uint16_t>::max ();

    const Expression &goal;
    const WalkNet<T> &net;
    // goal place -> hop table indexed by place
    std::unordered_map<size_t, std::vector<uint16_t>> hops;
    size_t maxMarkedForStructural; // beyond this many marked places, fall back

    void collectGoalPlaces (const Expression &e, std::vector<size_t> &out) const
    {
      if (e.kind == Expression::Kind::Atom) {
        const auto &a = e.atom;
        if (a.terms.size () == 1 && a.terms[0].second == 1 && a.constant >= 1
            && (a.op == petri::expr::Cmp::GE || a.op == petri::expr::Cmp::EQ)) {
          out.push_back (a.terms[0].first);
        }
        return;
      }
      for (const auto &c : e.children) collectGoalPlaces (c, out);
    }

    /** Backward BFS from g: hops[q] = transitions from a token in q to g. */
    std::vector<uint16_t> computeHops (size_t g, const MatrixCol<T> &producers) const
    {
      std::vector<uint16_t> d (net.placeCount (), UNREACHED);
      std::vector<size_t> queue;
      d[g] = 0;
      queue.push_back (g);
      for (size_t qi = 0; qi < queue.size (); ++qi) {
        size_t q = queue[qi];
        uint16_t nd = static_cast<uint16_t> (d[q] + 1);
        const SparseArray<T> &prods = producers.getColumn (q);
        for (size_t i = 0; i < prods.size (); ++i) {
          const SparseArray<T> &pre = net.pre (prods.keyAt (i));
          for (size_t j = 0; j < pre.size (); ++j) {
            size_t r = pre.keyAt (j);
            if (d[r] == UNREACHED) {
              d[r] = nd;
              queue.push_back (r);
            }
          }
        }
      }
      return d;
    }

    uint64_t atomDistance (const petri::expr::LinearAtom &a, const Marking<T> &m) const
    {
      long long v = a.value (m);
      uint64_t md = petri::expr::atomDistance (v, a.op, a.constant);
      if (md == 0 || a.terms.size () != 1) return md;
      auto it = hops.find (a.terms[0].first);
      if (it == hops.end ()) return md;
      const SparseArray<T> &sp = m.sparse ();
      if (sp.size () > maxMarkedForStructural) return md;
      const std::vector<uint16_t> &h = it->second;
      uint64_t best = petri::expr::INFINITE_DISTANCE;
      for (size_t i = 0; i < sp.size (); ++i) {
        uint16_t hq = h[sp.keyAt (i)];
        if (hq < best) best = hq;
      }
      // each missing token needs at least one journey; keep md as a floor
      return best >= petri::expr::INFINITE_DISTANCE ? md : (best > md ? best : md);
    }

    uint64_t distanceRec (const Expression &e, const Marking<T> &m) const
    {
      switch (e.kind) {
      case Expression::Kind::True: return 0;
      case Expression::Kind::False: return petri::expr::INFINITE_DISTANCE;
      case Expression::Kind::Not: return e.eval (m) ? 0 : 1;
      case Expression::Kind::And: {
        uint64_t sum = 0;
        for (const auto &c : e.children) {
          uint64_t d = distanceRec (c, m);
          if (d >= petri::expr::INFINITE_DISTANCE) return petri::expr::INFINITE_DISTANCE;
          sum += d;
        }
        return sum;
      }
      case Expression::Kind::Or: {
        uint64_t best = petri::expr::INFINITE_DISTANCE;
        for (const auto &c : e.children) {
          uint64_t d = distanceRec (c, m);
          if (d < best) best = d;
          if (best == 0) break;
        }
        return best;
      }
      case Expression::Kind::Atom: return atomDistance (e.atom, m);
      }
      return petri::expr::INFINITE_DISTANCE;
    }

  public:
    StructuralDistance (const Expression &g, const WalkNet<T> &n, size_t maxMarked = 256)
        : goal (g), net (n), maxMarkedForStructural (maxMarked)
    {
      std::vector<size_t> places;
      collectGoalPlaces (goal, places);
      if (places.empty ()) return;
      MatrixCol<T> producers = net.getNet ().getFlowTP ().transpose ();
      for (size_t p : places) {
        if (hops.find (p) == hops.end ()) hops.emplace (p, computeHops (p, producers));
      }
    }

    size_t goalPlaceCount () const
    {
      return hops.size ();
    }

    uint64_t of (const Marking<T> &m) const override
    {
      return distanceRec (goal, m);
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_STRUCTURALDISTANCE_H_ */
