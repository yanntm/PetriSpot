/*
 * DeadTransitions.h
 *
 * Transitions the state equation proves never enabled, places first. The
 * base program is m0 + C.x >= 0 (one row per place some live transition
 * changes), x >= 0. Phase one maximises each place's marking over it: a
 * bound below one means the place never holds a token, so every transition
 * touching it is dead (a producer would have marked it); a finite bound
 * kills the consumers asking more than it. Phase two, on what survives:
 * m0[p] + C[p].x >= pre(t)[p] for each input place of t as extra rows,
 * infeasible over the rationals means no reachable marking enables t. One
 * solve per place, then one per remaining transition, on the sparse edit
 * state of a reduction (flow matrices and their transposes, marking, live
 * flags), under a deadline and a pivot cap per solve. A net declared
 * one-safe adds m0[p] + C[p].x <= 1 to every place row: what the invariant
 * set of libHSC's linear test knows, the state equation now knows too.
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
  size_t places = 0, stuck = 0, byBound = 0; // phase one: places bounded, places never marked, transitions dead by a bound
  long placesMs = 0;                          // phase one's time
  bool limited = false; // the deadline stopped the sweep
};

template<typename T>
  std::vector<size_t> deadTransitions (const MatrixCol<T> &pre, const MatrixCol<T> &post,
                                       const MatrixCol<T> &consumers, const MatrixCol<T> &producers,
                                       const std::vector<T> &marks, const std::vector<bool> &liveT,
                                       const std::vector<bool> &liveP, bool safe,
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
      if (safe) r.hi = 1.0 - static_cast<double> (marks[p]);
      rowOf[p] = base.addRow (std::move (r));
    }
    for (size_t t = 0; t < nt; ++t) if (!liveT[t]) base.upper[t] = 0.0; // a retired slot never fires
    LpLimits limits;
    limits.hasDeadline = true;
    limits.deadline = deadline;
    limits.maxPivots = maxPivots;
    const auto phaseOne = std::chrono::steady_clock::now ();
    std::vector<bool> gone (nt, false); // dead by phase one, out of phase two and fixed at zero in its programs
    auto kill = [&] (size_t t) {
      if (gone[t] || !liveT[t]) return;
      gone[t] = true;
      dead.push_back (t);
      base.upper[t] = 0.0;
      ++report.dead;
      ++report.byBound;
    };
    for (size_t p = 0; p < np; ++p) {
      if (!liveP[p] || rowOf[p] == absent) continue;
      if (std::chrono::steady_clock::now () >= deadline) { report.limited = true; return dead; }
      std::vector<double> c (nt, 0.0);
      const SparseArray<long long> &eff = byPlace.getColumn (p);
      for (size_t i = 0; i < eff.size (); ++i) c[eff.keyAt (i)] = static_cast<double> (eff.valueAt (i));
      base.setMaximise (c);
      Simplex simplex (limits);
      LpResult r = simplex.solve (base);
      ++report.solves;
      report.pivots += r.pivots;
      ++report.places;
      if (r.status == LpStatus::TimeLimit) { report.limited = true; base.objective.clear (); return dead; }
      if (r.status != LpStatus::Optimal) continue; // unbounded, or the pivot cap: no bound
      const double bound = static_cast<double> (marks[p]) - r.objective; // the objective is the negated maximum
      if (bound < 1.0 - 1e-6) {
        ++report.stuck;
        for (const auto *mat : { &consumers, &producers }) {
          const SparseArray<T> &ts = mat->getColumn (p);
          for (size_t i = 0; i < ts.size (); ++i) kill (ts.keyAt (i));
        }
      } else {
        const SparseArray<T> &ts = consumers.getColumn (p);
        for (size_t i = 0; i < ts.size (); ++i)
          if (static_cast<double> (ts.valueAt (i)) > bound + 1e-6) kill (ts.keyAt (i));
      }
    }
    base.objective.clear ();
    report.placesMs = std::chrono::duration_cast<std::chrono::milliseconds> (std::chrono::steady_clock::now () - phaseOne).count ();
    const size_t start = cursor < nt ? cursor : 0;
    for (size_t k = 0; k < nt; ++k) {
      const size_t t = (start + k) % nt;
      if (!liveT[t] || gone[t]) continue;
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
