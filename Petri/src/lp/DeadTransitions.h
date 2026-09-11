/*
 * DeadTransitions.h
 *
 * Transitions the state equation proves never enabled. For a transition t
 * the program is m0 + C.x >= 0 (one row per place some live transition
 * changes), x >= 0, and m0[p] + C[p].x >= pre(t)[p] for each input place of
 * t; infeasible over the rationals means no reachable marking enables t.
 * Solved on the sparse edit state of a reduction (flow matrices, marking,
 * live flags), one base problem and one solve with the enabling rows as
 * extra rows per transition, under a deadline and a pivot cap per solve.
 * A verdict is the floating-point solver's; the exact Farkas checker of
 * algorithm.md is still to come. `cursor` is where the sweep starts and
 * where it stopped, so a pass cut by its budget resumes after the last
 * transition tested.
 */
#ifndef PETRI_LP_DEADTRANSITIONS_H_
#define PETRI_LP_DEADTRANSITIONS_H_
#include <chrono>
#include <limits>
#include <vector>
#include "core/MatrixCol.h"
#include "core/SparseArray.h"
#include "lp/LpProblem.h"
#include "lp/Simplex.h"

namespace petri::lp
{

struct DeadReport
{
  size_t tested = 0, dead = 0, alive = 0, untested = 0, solves = 0, pivots = 0;
  bool limited = false; // the deadline stopped the sweep
};

template<typename T>
  std::vector<size_t> deadTransitions (const MatrixCol<T> &pre, const MatrixCol<T> &post,
                                       const std::vector<T> &marks, const std::vector<bool> &liveT,
                                       size_t &cursor, std::chrono::steady_clock::time_point deadline,
                                       size_t maxPivots, DeadReport &report)
  {
    constexpr size_t absent = std::numeric_limits<size_t>::max ();
    const size_t nt = pre.getColumnCount (), np = marks.size ();
    std::vector<size_t> dead;
    if (nt == 0) return dead;
    // the effect of every live transition on each place: the rows of the program
    MatrixCol<long long> byPlace (nt, np);
    for (size_t t = 0; t < nt; ++t) {
      if (!liveT[t]) continue;
      const SparseArray<T> &in = pre.getColumn (t);
      const SparseArray<T> &out = post.getColumn (t);
      for (size_t i = 0; i < in.size (); ++i) {
        SparseArray<long long> &col = byPlace.getColumn (in.keyAt (i));
        col.put (t, col.get (t) - static_cast<long long> (in.valueAt (i)));
      }
      for (size_t i = 0; i < out.size (); ++i) {
        SparseArray<long long> &col = byPlace.getColumn (out.keyAt (i));
        col.put (t, col.get (t) + static_cast<long long> (out.valueAt (i)));
      }
    }
    LpProblem base (nt);
    std::vector<size_t> rowOf (np, absent);
    for (size_t p = 0; p < np; ++p) {
      const SparseArray<long long> &col = byPlace.getColumn (p);
      if (col.size () == 0) continue;
      Row r;
      r.coeffs = col;
      r.lo = -static_cast<double> (marks[p]);
      rowOf[p] = base.addRow (std::move (r));
    }
    for (size_t t = 0; t < nt; ++t) if (!liveT[t]) base.upper[t] = 0.0; // a retired slot never fires
    LpLimits limits;
    limits.hasDeadline = true;
    limits.deadline = deadline;
    limits.maxPivots = maxPivots;
    const size_t start = cursor < nt ? cursor : 0;
    for (size_t k = 0; k < nt; ++k) {
      const size_t t = (start + k) % nt;
      if (!liveT[t]) continue;
      if (std::chrono::steady_clock::now () >= deadline) { report.limited = true; cursor = t; return dead; }
      const SparseArray<T> &in = pre.getColumn (t);
      std::vector<Row> extra;
      bool constantDead = false;
      for (size_t i = 0; i < in.size (); ++i) {
        const size_t p = in.keyAt (i);
        const T w = in.valueAt (i);
        if (rowOf[p] == absent) { if (marks[p] < w) constantDead = true; continue; } // a place nothing changes
        Row r;
        r.coeffs = byPlace.getColumn (p);
        r.lo = static_cast<double> (w) - static_cast<double> (marks[p]);
        extra.push_back (std::move (r));
      }
      ++report.tested;
      if (constantDead) { dead.push_back (t); ++report.dead; continue; }
      if (extra.empty ()) { ++report.alive; continue; } // enabled at m0, or no input place
      Simplex simplex (limits);
      LpResult r = simplex.solve (base, extra);
      ++report.solves;
      report.pivots += r.pivots;
      if (r.status == LpStatus::Infeasible) { dead.push_back (t); ++report.dead; }
      else if (r.feasible ()) ++report.alive;
      else if (r.status == LpStatus::TimeLimit) { report.limited = true; --report.tested; cursor = t; return dead; }
      else ++report.untested; // pivot cap or too large: this one stays, the sweep goes on
    }
    cursor = start; // a full sweep: the next starts where this one did
    return dead;
  }

} // namespace petri::lp
#endif /* PETRI_LP_DEADTRANSITIONS_H_ */
