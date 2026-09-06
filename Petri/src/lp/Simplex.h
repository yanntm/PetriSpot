/*
 * Simplex.h
 *
 * A bounded-variable primal simplex in double precision over the sparse
 * columns of an LpProblem. Every row i becomes the equality
 * a_i.x - r_i + art_i = 0 with the slack r_i bounded by [lo_i, hi_i] and an
 * artificial art_i that phase 1 drives to zero (minimising the sum of their
 * absolute values); phase 2 minimises the objective with the artificials
 * fixed at zero. The basis inverse is kept in product form (Basis.h), one
 * sparse eta per pivot, and rebuilt from the basis columns periodically.
 * Pricing is Dantzig over a candidate list refreshed by full sweeps, the ratio test is
 * Harris's two-pass, Bland's rule takes over after a run of degenerate
 * pivots. Limits: a pivot budget, a deadline, a row count the dense inverse
 * can afford. Knows nothing of Petri nets. See algorithm.md section 3.
 */
#ifndef PETRI_LP_SIMPLEX_H_
#define PETRI_LP_SIMPLEX_H_

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

#include "core/MatrixCol.h"
#include "core/SparseArray.h"
#include "lp/Basis.h"
#include "lp/LpProblem.h"

namespace petri::lp
{

enum class LpStatus
{
  Optimal, Infeasible, Unbounded, PivotLimit, TimeLimit, TooLarge
};

inline const char* to_string (LpStatus s)
{
  switch (s) {
  case LpStatus::Optimal: return "optimal";
  case LpStatus::Infeasible: return "infeasible";
  case LpStatus::Unbounded: return "unbounded";
  case LpStatus::PivotLimit: return "pivot limit";
  case LpStatus::TimeLimit: return "time limit";
  case LpStatus::TooLarge: return "too large";
  }
  return "?";
}

struct LpLimits
{
  size_t maxPivots = 2000000;
  std::chrono::steady_clock::time_point deadline;
  bool hasDeadline = false;
  size_t maxRows = std::numeric_limits<size_t>::max (); // no limit: the inverse is sparse
};

struct LpResult
{
  LpStatus status = LpStatus::TooLarge;
  std::vector<double> x;   // structural columns, when a solution exists
  double objective = 0.0;  // of the problem's objective, when optimal with one
  size_t pivots = 0;
  size_t infeasibleRow = std::numeric_limits<size_t>::max (); // a row whose artificial stayed positive
  bool feasible () const
  {
    return status == LpStatus::Optimal || status == LpStatus::Unbounded;
  }
};

class Simplex
{
  static constexpr double TOL_FEAS = 1e-7;
  static constexpr double TOL_OPT = 1e-9;
  static constexpr double TOL_PIVOT = 1e-9;
  static constexpr size_t REFACTOR = 200;     // pivots between two refactorisations
  static constexpr size_t DEGENERATE = 64;    // zero-step pivots before Bland's rule
  static constexpr uint32_t NONE = std::numeric_limits<uint32_t>::max ();

  size_t m = 0, nStruct = 0, n = 0;          // rows, structural columns, all columns
  const MatrixCol<long long> *cols = nullptr; // structural columns of the base problem, by row
  MatrixCol<long long> extraCols;             // the extra rows of this solve, by column; their row index is mBase + k
  size_t mBase = 0;                           // rows of the base problem
  std::vector<double> lower, upper, cost, x;  // per column (structural, slacks, artificials)
  std::vector<uint32_t> basis;                // per row, the basic column
  std::vector<int32_t> rowOf;                 // per column, its row when basic, -1 otherwise
  Basis binv;                                 // the basis inverse in product form
  std::vector<double> work;                   // a dense vector of size m for ftran
  std::vector<double> y, alpha;               // duals of the phase, entering column in the basis
  std::vector<uint32_t> candidates;           // multiple pricing: the best columns of the last sweep
  static constexpr size_t CANDIDATES = 64;
  size_t degenerate = 0;
  bool bland = false;
  size_t pivots = 0, sinceRefactor = 0;
  double minPivot = INF; size_t blandPivots = 0; // debug statistics of the pivots since the last refactor
  LpLimits limits;

public:
  bool debug = false;
  bool reported = false;

private:

  template<class F>
    void forEachEntry (uint32_t j, F f) const
    {
      if (j < nStruct) {
        const SparseArray<long long> &c = cols->getColumn (j);
        for (size_t i = 0, e = c.size (); i < e; ++i) f (static_cast<uint32_t> (c.keyAt (i)), static_cast<double> (c.valueAt (i)));
        if (mBase < m) {
          const SparseArray<long long> &x2 = extraCols.getColumn (j);
          for (size_t i = 0, e = x2.size (); i < e; ++i) f (static_cast<uint32_t> (mBase + x2.keyAt (i)), static_cast<double> (x2.valueAt (i)));
        }
      } else if (j < nStruct + m) {
        f (static_cast<uint32_t> (j - nStruct), -1.0);
      } else {
        f (static_cast<uint32_t> (j - nStruct - m), 1.0);
      }
    }

  bool atLower (uint32_t j) const
  {
    return lower[j] > -INF && x[j] <= lower[j] + TOL_FEAS;
  }
  bool atUpper (uint32_t j) const
  {
    return upper[j] < INF && x[j] >= upper[j] - TOL_FEAS;
  }

  /**
   * The eta file from scratch. A basic unit column (a slack -e_i or an
   * artificial +e_i, wherever the pivots left it) is assigned to its own
   * row i and goes into the signed diagonal; a set of basic columns has no
   * row order of its own. The remaining rows are free, and each structural
   * basic column is pivoted into whichever free row gives it the largest
   * pivot, sparsest columns first. A column that finds no usable pivot made
   * the basis singular: it leaves the basis at its current value and its row
   * takes the artificial. Then the basic values are recomputed from the
   * nonbasic ones.
   */
  void refactor ()
  {
    std::vector<double> d (m, 1.0);
    std::vector<uint32_t> newBasis (m, NONE);
    std::vector<uint32_t> structural;
    for (size_t k = 0; k < m; ++k) {
      uint32_t b = basis[k];
      if (b < nStruct) { structural.push_back (b); continue; }
      size_t i = b < nStruct + m ? b - nStruct : b - nStruct - m;
      if (newBasis[i] != NONE) { structural.push_back (b); continue; } // the slack and the artificial of one row: dependent
      newBasis[i] = b;
      d[i] = b < nStruct + m ? -1.0 : 1.0;
    }
    std::vector<size_t> freeRows;
    for (size_t r = 0; r < m; ++r)
      if (newBasis[r] == NONE) freeRows.push_back (r);
    binv.reset (m, d);
    std::sort (structural.begin (), structural.end (), [&] (uint32_t a, uint32_t b) {
      return cols->getColumn (a).size () < cols->getColumn (b).size ();
    });
    size_t dropped = 0;
    for (uint32_t j : structural) {
      std::fill (work.begin (), work.end (), 0.0);
      forEachEntry (j, [&] (uint32_t r, double v) { work[r] = v; });
      binv.ftran (work);
      size_t bestRow = m;
      double bestPivot = TOL_PIVOT;
      for (size_t r : freeRows)
        if (std::fabs (work[r]) > bestPivot) { bestPivot = std::fabs (work[r]); bestRow = r; }
      if (bestRow == m) { ++dropped; continue; } // dependent: out of the basis, its row keeps the artificial
      binv.pivot (bestRow, work);
      newBasis[bestRow] = j;
      freeRows.erase (std::find (freeRows.begin (), freeRows.end (), bestRow));
    }
    for (size_t r : freeRows) newBasis[r] = static_cast<uint32_t> (nStruct + m + r);
    basis = std::move (newBasis);
    std::fill (rowOf.begin (), rowOf.end (), -1);
    for (size_t k = 0; k < m; ++k) rowOf[basis[k]] = static_cast<int32_t> (k);
    // x_B = -B^-1 (sum over nonbasic j of A_j x_j)
    std::fill (work.begin (), work.end (), 0.0);
    for (uint32_t j = 0; j < n; ++j) {
      if (rowOf[j] >= 0 || x[j] == 0.0) continue;
      forEachEntry (j, [&] (uint32_t i, double v) { work[i] += v * x[j]; });
    }
    binv.ftran (work);
    for (size_t k = 0; k < m; ++k) x[basis[k]] = -work[k];
    sinceRefactor = 0;
    if (debug) std::cerr << "refactor at pivot " << pivots << ": " << structural.size () - dropped << " structural basics, "
        << dropped << " dropped, etas " << binv.etaCount () << " (" << binv.nonZeros () << " entries), residual " << residual ()
        << "; since the last: smallest pivot " << minPivot << ", " << blandPivots << " Bland pivots" << std::endl;
    minPivot = INF;
    blandPivots = 0;
  }

  /** Duals of the current cost: y = c_B B^-1. */
  void computeDuals (const std::vector<double> &c)
  {
    for (size_t k = 0; k < m; ++k) y[k] = c[basis[k]];
    binv.btran (y);
  }

  double reducedCost (uint32_t j, const std::vector<double> &c) const
  {
    double d = c[j];
    forEachEntry (j, [&] (uint32_t i, double v) { d -= y[i] * v; });
    return d;
  }

  /** Whether column j may enter, and in which direction (+1 increases, -1 decreases); 0 otherwise. */
  int eligible (uint32_t j, const std::vector<double> &c, double &d) const
  {
    if (rowOf[j] >= 0 || lower[j] == upper[j]) return 0;
    d = reducedCost (j, c);
    if (d < -TOL_OPT && !atUpper (j)) return 1;
    if (d > TOL_OPT && !atLower (j)) return -1;
    return 0;
  }

  /**
   * The entering column and its direction; NONE when the basis is optimal.
   * Multiple pricing: a full sweep keeps the CANDIDATES columns of largest
   * reduced cost, and the following iterations choose among them (Dantzig)
   * until none is eligible, when the sweep runs again. Under Bland's rule the
   * columns are scanned from the first and the first eligible one enters.
   */
  uint32_t price (const std::vector<double> &c, size_t priced, int &dir)
  {
    double d;
    if (bland) {
      for (uint32_t col = 0; col < priced; ++col) {
        int dj = eligible (col, c, d);
        if (dj != 0) { dir = dj; return col; }
      }
      return NONE;
    }
    uint32_t best = NONE;
    double bestScore = 0.0;
    int bestDir = 0;
    for (uint32_t col : candidates) {
      int dj = eligible (col, c, d);
      if (dj != 0 && std::fabs (d) > bestScore) { best = col; bestScore = std::fabs (d); bestDir = dj; }
    }
    if (best != NONE) { dir = bestDir; return best; }
    candidates.clear ();
    std::vector<std::pair<double, uint32_t>> found; // (-|d|, column)
    for (uint32_t col = 0; col < priced; ++col) {
      int dj = eligible (col, c, d);
      if (dj == 0) continue;
      found.emplace_back (-std::fabs (d), col);
      if (std::fabs (d) > bestScore) { best = col; bestScore = std::fabs (d); bestDir = dj; }
    }
    if (!found.empty ()) {
      size_t keep = std::min (found.size (), CANDIDATES);
      std::partial_sort (found.begin (), found.begin () + static_cast<long> (keep), found.end ());
      for (size_t i = 0; i < keep; ++i) candidates.push_back (found[i].second);
    }
    dir = bestDir;
    return best;
  }

  /** alpha = B^-1 A_q. */
  void enteringColumn (uint32_t q)
  {
    std::fill (alpha.begin (), alpha.end (), 0.0);
    forEachEntry (q, [&] (uint32_t i, double v) { alpha[i] = v; });
    binv.ftran (alpha);
  }

  /** Harris two-pass ratio test; returns the leaving row, m for a bound flip, NONE-like m+1 when unbounded. */
  size_t ratioTest (uint32_t q, int dir, double &theta)
  {
    double flip = (lower[q] > -INF && upper[q] < INF) ? upper[q] - lower[q] : INF;
    double thetaMax = flip;
    for (size_t k = 0; k < m; ++k) {
      double a = alpha[k];
      if (std::fabs (a) < TOL_PIVOT) continue;
      uint32_t b = basis[k];
      double delta = -dir * a; // change of x_b per unit of theta
      double t;
      if (delta < 0 && lower[b] > -INF) t = (x[b] - lower[b] + TOL_FEAS) / -delta;
      else if (delta > 0 && upper[b] < INF) t = (upper[b] - x[b] + TOL_FEAS) / delta;
      else continue;
      if (t < thetaMax) thetaMax = t;
    }
    if (thetaMax == INF) return m + 1;
    size_t leave = m;
    double bestPivot = 0.0;
    theta = flip;
    for (size_t k = 0; k < m; ++k) {
      double a = alpha[k];
      if (std::fabs (a) < TOL_PIVOT) continue;
      uint32_t b = basis[k];
      double delta = -dir * a;
      double t;
      if (delta < 0 && lower[b] > -INF) t = (x[b] - lower[b]) / -delta;
      else if (delta > 0 && upper[b] < INF) t = (upper[b] - x[b]) / delta;
      else continue;
      if (t > thetaMax) continue;
      // pass 2: among the rows within the relaxed step, the largest pivot; under Bland, the smallest basic index
      bool better = leave == m || (bland ? b < basis[leave] : std::fabs (a) > bestPivot);
      if (better) {
        leave = k;
        bestPivot = std::fabs (a);
        theta = std::max (0.0, t);
      }
    }
    if (leave == m) theta = flip;
    return leave;
  }

  /** Move by theta along the entering direction, then pivot or flip. */
  void update (uint32_t q, int dir, double theta, size_t leave)
  {
    for (size_t k = 0; k < m; ++k)
      if (alpha[k] != 0.0) x[basis[k]] -= dir * theta * alpha[k];
    x[q] += dir * theta;
    if (leave == m) {
      x[q] = dir > 0 ? upper[q] : lower[q]; // a bound flip, the basis stands
      return;
    }
    uint32_t b = basis[leave];
    double delta = -dir * alpha[leave];
    x[b] = delta < 0 ? lower[b] : upper[b];
    if (debug) {
      minPivot = std::min (minPivot, std::fabs (alpha[leave]));
      if (bland) ++blandPivots;
    }
    binv.pivot (leave, alpha);
    if (debug && !reported) {
      // the new inverse must send the entering column to e_leave
      std::fill (work.begin (), work.end (), 0.0);
      forEachEntry (q, [&] (uint32_t r, double v) { work[r] = v; });
      binv.ftran (work);
      double dev = 0.0;
      for (size_t r = 0; r < m; ++r) dev = std::max (dev, std::fabs (work[r] - (r == leave ? 1.0 : 0.0)));
      if (dev > 1e-6) {
        reported = true;
        std::cerr << "pivot " << pivots << ": entering " << q << (q < nStruct ? " structural" : (q < nStruct + m ? " slack" : " artificial"))
            << " leaves row " << leave << " (was " << b << "), alpha[leave] " << alpha[leave] << ", theta " << theta
            << ", bland " << bland << ", etas " << binv.etaCount () << ", deviation " << dev << std::endl;
      }
    }
    rowOf[b] = -1;
    rowOf[q] = static_cast<int32_t> (leave);
    basis[leave] = q;
    ++pivots;
    ++sinceRefactor;
    if (sinceRefactor >= REFACTOR) {
      if (debug) {
        // does the maintained inverse still send every basic column to its unit vector?
        double worst = 0.0;
        size_t worstRow = m;
        for (size_t k = 0; k < m; ++k) {
          std::fill (work.begin (), work.end (), 0.0);
          forEachEntry (basis[k], [&] (uint32_t r, double v) { work[r] = v; });
          binv.ftran (work);
          for (size_t r = 0; r < m; ++r) {
            double dev = std::fabs (work[r] - (r == k ? 1.0 : 0.0));
            if (dev > worst) { worst = dev; worstRow = k; }
          }
        }
        std::cerr << "before refactor at pivot " << pivots << ": inverse deviation " << worst << " at row " << worstRow
            << " (basic column " << (worstRow < m ? basis[worstRow] : 0) << "), residual " << residual () << std::endl;
      }
      refactor ();
    }
  }

  /** Debug: the largest violation of a row equation a.x - r + art = 0 by the current point. */
  double residual () const
  {
    std::vector<double> v (m, 0.0);
    for (uint32_t j = 0; j < n; ++j) {
      if (x[j] == 0.0) continue;
      forEachEntry (j, [&] (uint32_t i, double a) { v[i] += a * x[j]; });
    }
    double worst = 0.0;
    for (size_t i = 0; i < m; ++i) worst = std::max (worst, std::fabs (v[i]));
    return worst;
  }

  bool outOfBudget (LpStatus &why) const
  {
    if (pivots >= limits.maxPivots) { why = LpStatus::PivotLimit; return true; }
    if (limits.hasDeadline && (pivots & 15) == 0 && std::chrono::steady_clock::now () >= limits.deadline) {
      why = LpStatus::TimeLimit;
      return true;
    }
    return false;
  }

  /** Iterate to optimality of cost c over the first `priced` columns. */
  LpStatus iterate (const std::vector<double> &c, size_t priced)
  {
    LpStatus why;
    for (;;) {
      if (outOfBudget (why)) return why;
      computeDuals (c);
      int dir = 0;
      uint32_t q = price (c, priced, dir);
      if (q == NONE) return LpStatus::Optimal;
      enteringColumn (q);
      double theta = 0.0;
      size_t leave = ratioTest (q, dir, theta);
      if (leave == m + 1) return LpStatus::Unbounded;
      if (theta <= TOL_FEAS) {
        if (++degenerate >= DEGENERATE) bland = true;
      } else {
        degenerate = 0;
        bland = false;
      }
      update (q, dir, theta, leave);
    }
  }

public:
  explicit Simplex (const LpLimits &lim = LpLimits ())
      : limits (lim)
  {
  }

  LpResult solve (const LpProblem &p)
  {
    return solve (p, std::vector<Row> ());
  }

  /** The base problem with extra rows appended for this solve only: a branch or a cut, without copying the base. */
  LpResult solve (const LpProblem &p, const std::vector<Row> &extra)
  {
    LpResult res;
    mBase = p.rowCount ();
    m = mBase + extra.size ();
    nStruct = p.columns;
    n = nStruct + 2 * m;
    if (m > limits.maxRows) {
      res.status = LpStatus::TooLarge;
      return res;
    }
    cols = &p.byColumn ();
    if (!extra.empty ()) {
      MatrixCol<long long> byRow (nStruct, 0);
      for (const Row &r : extra) byRow.appendColumn (r.coeffs);
      extraCols = byRow.transpose ();
    }
    lower.assign (n, 0.0);
    upper.assign (n, INF);
    cost.assign (n, 0.0);
    x.assign (n, 0.0);
    for (size_t j = 0; j < nStruct; ++j) {
      lower[j] = p.lower[j];
      upper[j] = p.upper[j];
      if (p.hasObjective ()) cost[j] = p.objective[j];
      x[j] = lower[j] > -INF ? lower[j] : (upper[j] < INF ? upper[j] : 0.0);
    }
    for (size_t i = 0; i < m; ++i) {
      lower[nStruct + i] = i < mBase ? p.lo[i] : extra[i - mBase].lo;
      upper[nStruct + i] = i < mBase ? p.hi[i] : extra[i - mBase].hi;
    }
    // row activities at the starting point, slacks clamped into their range, artificials take the rest
    std::vector<double> v (m, 0.0);
    for (size_t j = 0; j < nStruct; ++j)
      if (x[j] != 0.0)
        forEachEntry (static_cast<uint32_t> (j), [&] (uint32_t i, double a) { v[i] += a * x[j]; });
    std::vector<double> cost1 (n, 0.0);
    basis.assign (m, 0);
    rowOf.assign (n, -1);
    for (size_t i = 0; i < m; ++i) {
      uint32_t s = static_cast<uint32_t> (nStruct + i), a = static_cast<uint32_t> (nStruct + m + i);
      x[s] = std::min (std::max (v[i], lower[s]), upper[s]);
      x[a] = x[s] - v[i];
      // a row satisfied at the start puts its slack in the basis, its artificial fixed at zero and costless;
      // a violated row starts on its artificial, which phase 1 drives out
      if (x[a] == 0.0) {
        lower[a] = upper[a] = 0.0;
        basis[i] = s;
        rowOf[s] = static_cast<int32_t> (i);
      } else {
        if (x[a] > 0) { lower[a] = 0; upper[a] = INF; cost1[a] = 1.0; }
        else { lower[a] = -INF; upper[a] = 0; cost1[a] = -1.0; }
        basis[i] = a;
        rowOf[a] = static_cast<int32_t> (i);
      }
    }
    // the basis matrix is diagonal (+1 artificial, -1 slack): its inverse is itself
    std::vector<double> d (m);
    for (size_t i = 0; i < m; ++i) d[i] = basis[i] < nStruct + m ? -1.0 : 1.0;
    binv.reset (m, d);
    y.assign (m, 0.0);
    alpha.assign (m, 0.0);
    work.assign (m, 0.0);
    candidates.clear ();
    degenerate = 0;
    bland = false;
    pivots = 0;
    sinceRefactor = 0;
    // phase 1: the artificials to zero
    LpStatus st = iterate (cost1, n);
    res.pivots = pivots;
    if (st != LpStatus::Optimal && st != LpStatus::Unbounded) {
      res.status = st;
      return res;
    }
    double infeas = 0.0;
    size_t worst = std::numeric_limits<size_t>::max ();
    double worstVal = 0.0;
    for (size_t i = 0; i < m; ++i) {
      double a = std::fabs (x[nStruct + m + i]);
      infeas += a;
      if (a > worstVal) { worstVal = a; worst = i; }
    }
    if (infeas > TOL_FEAS * static_cast<double> (m + 1)) {
      res.status = LpStatus::Infeasible;
      res.infeasibleRow = worst;
      return res;
    }
    // the artificials are fixed at zero from here on
    for (size_t i = 0; i < m; ++i) {
      uint32_t a = static_cast<uint32_t> (nStruct + m + i);
      lower[a] = upper[a] = 0.0;
      if (rowOf[a] < 0) x[a] = 0.0;
    }
    if (p.hasObjective ()) {
      degenerate = 0;
      bland = false;
      candidates.clear ();
      st = iterate (cost, nStruct + m);
      res.pivots = pivots;
      if (st != LpStatus::Optimal && st != LpStatus::Unbounded) {
        res.status = st;
        res.x.assign (x.begin (), x.begin () + static_cast<long> (nStruct)); // feasible, not optimal
        return res;
      }
      res.status = st;
    } else {
      res.status = LpStatus::Optimal;
    }
    res.x.assign (x.begin (), x.begin () + static_cast<long> (nStruct));
    res.objective = 0.0;
    if (p.hasObjective ())
      for (size_t j = 0; j < nStruct; ++j) res.objective += p.objective[j] * res.x[j];
    return res;
  }
};

} // namespace petri::lp

#endif /* PETRI_LP_SIMPLEX_H_ */
