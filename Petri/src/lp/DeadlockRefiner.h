/*
 * DeadlockRefiner.h
 *
 * A dead marking has, for every transition, a pre-place below the weight of
 * the arc: a conjunction over the transitions of one disjunction each, far
 * too many conjuncts for a normal form. This refiner builds it lazily: it
 * reads the marking of a candidate solution, finds a transition still
 * enabled there, and splits the problem on which of its pre-places to
 * starve (one branch per pre-place, the row s_p <= w - 1). Every row added
 * is a necessary condition of a dead marking, so a problem whose every
 * branch is infeasible has no reachable deadlock, up to the relaxation.
 * The transition with the fewest pre-places is split first, the pre-place
 * with the least tokens in the candidate first. See algorithm.md section 4.
 */
#ifndef PETRI_LP_DEADLOCKREFINER_H_
#define PETRI_LP_DEADLOCKREFINER_H_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

#include "core/MatrixCol.h"
#include "core/SparseArray.h"
#include "core/SparsePetriNet.h"
#include "lp/LpProblem.h"
#include "lp/Refiner.h"
#include "lp/Simplex.h"
#include "lp/StateEquation.h"

namespace petri::lp
{

template<typename T>
  class DeadlockRefiner : public Refiner
  {
    const SparsePetriNet<T> &net;
    const MatrixCol<T> &effects;  // by transition
    MatrixCol<T> effectsT;        // by place
    std::vector<double> marking;  // of the candidate, per place
    static constexpr double TOL = 1e-6;

  public:
    size_t splits = 0;

    DeadlockRefiner (const SparsePetriNet<T> &n, const MatrixCol<T> &eff)
        : net (n), effects (eff), effectsT (eff.transpose ())
    {
    }

    /** The row s_p <= w - 1 over the transition counts: sum_t C[p,t] x_t <= w - 1 - m0_p. */
    Row starve (size_t p, T w) const
    {
      Row r;
      r.coeffs = StateEquation<T>::toLongLong (effectsT.getColumn (p));
      r.lo = -INF;
      r.hi = static_cast<double> (w) - 1.0 - static_cast<double> (net.getMarks ()[p]);
      return r;
    }

    Verdict examine (const LpProblem &, const LpResult &solution) override
    {
      // the candidate marking
      marking.assign (net.getPlaceCount (), 0.0);
      for (size_t p = 0; p < marking.size (); ++p) marking[p] = static_cast<double> (net.getMarks ()[p]);
      for (size_t t = 0; t < solution.x.size (); ++t) {
        double k = solution.x[t];
        if (k <= TOL) continue;
        const SparseArray<T> &col = effects.getColumn (t);
        for (size_t i = 0; i < col.size (); ++i) marking[col.keyAt (i)] += k * static_cast<double> (col.valueAt (i));
      }
      // an enabled transition with the fewest pre-places
      size_t best = std::numeric_limits<size_t>::max (), bestArity = std::numeric_limits<size_t>::max ();
      for (size_t t = 0; t < net.getTransitionCount (); ++t) {
        const SparseArray<T> &pre = net.getFlowPT ().getColumn (t);
        if (pre.size () >= bestArity) continue;
        bool enabled = true;
        for (size_t i = 0; i < pre.size () && enabled; ++i)
          if (marking[pre.keyAt (i)] < static_cast<double> (pre.valueAt (i)) - TOL) enabled = false;
        if (enabled) {
          best = t;
          bestArity = pre.size ();
          if (bestArity <= 1) break;
        }
      }
      Verdict v;
      if (best == std::numeric_limits<size_t>::max ()) return v; // dead: accept the candidate
      const SparseArray<T> &pre = net.getFlowPT ().getColumn (best);
      if (pre.size () == 0) {
        // a transition without pre-places is always enabled: no dead marking exists; an empty split is that
        v.kind = Verdict::Kind::Split;
        return v;
      }
      std::vector<std::pair<double, size_t>> order; // (tokens, index in pre): the emptiest place first
      for (size_t i = 0; i < pre.size (); ++i) order.emplace_back (marking[pre.keyAt (i)], i);
      std::sort (order.begin (), order.end ());
      v.kind = Verdict::Kind::Split;
      for (const auto &o : order) v.branches.push_back ({ starve (pre.keyAt (o.second), pre.valueAt (o.second)) });
      // the loop pops the last branch first: push in reverse so the emptiest place is explored first
      std::reverse (v.branches.begin (), v.branches.end ());
      ++splits;
      return v;
    }
  };

} // namespace petri::lp

#endif /* PETRI_LP_DEADLOCKREFINER_H_ */
