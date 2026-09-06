/*
 * Simplex.h
 *
 * A bounded-variable primal simplex in double precision over the sparse
 * columns of an LpProblem. Every row i becomes the equality
 * a_i.x - r_i + art_i = 0 with the slack r_i bounded by [lo_i, hi_i] and an
 * artificial art_i that phase 1 drives to zero (minimising the sum of their
 * absolute values); phase 2 minimises the objective with the artificials
 * fixed at zero. The basis inverse is dense, updated by the product form at
 * each pivot and refactorised periodically by Gauss-Jordan elimination.
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
#include <limits>
#include <utility>
#include <vector>

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
  size_t maxRows = 6000; // the dense inverse: maxRows^2 doubles
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
  std::vector<std::vector<std::pair<uint32_t, double>>> cols; // structural columns by row
  std::vector<double> lower, upper, cost, x;  // per column (structural, slacks, artificials)
  std::vector<uint32_t> basis;                // per row, the basic column
  std::vector<int32_t> rowOf;                 // per column, its row when basic, -1 otherwise
  std::vector<double> Binv;                   // m x m, row-major
  std::vector<double> y, alpha;               // duals of the phase, entering column in the basis
  std::vector<uint32_t> candidates;           // multiple pricing: the best columns of the last sweep
  static constexpr size_t CANDIDATES = 64;
  size_t degenerate = 0;
  bool bland = false;
  size_t pivots = 0, sinceRefactor = 0;
  LpLimits limits;

  template<class F>
    void forEachEntry (uint32_t j, F f) const
    {
      if (j < nStruct) {
        for (const auto &e : cols[j]) f (e.first, e.second);
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

  /** Binv from scratch: Gauss-Jordan on the basis columns; then the basic values from the nonbasic ones. */
  void refactor ()
  {
    std::vector<double> B (m * m, 0.0);
    for (size_t k = 0; k < m; ++k)
      forEachEntry (basis[k], [&] (uint32_t i, double v) { B[i * m + k] = v; });
    std::fill (Binv.begin (), Binv.end (), 0.0);
    for (size_t i = 0; i < m; ++i) Binv[i * m + i] = 1.0;
    for (size_t c = 0; c < m; ++c) {
      size_t p = c;
      for (size_t r = c + 1; r < m; ++r)
        if (std::fabs (B[r * m + c]) > std::fabs (B[p * m + c])) p = r;
      if (std::fabs (B[p * m + c]) < 1e-12) continue; // singular: keep going, the next pivot repairs
      if (p != c) {
        for (size_t k = 0; k < m; ++k) {
          std::swap (B[p * m + k], B[c * m + k]);
          std::swap (Binv[p * m + k], Binv[c * m + k]);
        }
      }
      double d = 1.0 / B[c * m + c];
      for (size_t k = 0; k < m; ++k) {
        B[c * m + k] *= d;
        Binv[c * m + k] *= d;
      }
      for (size_t r = 0; r < m; ++r) {
        if (r == c) continue;
        double f = B[r * m + c];
        if (f == 0.0) continue;
        for (size_t k = 0; k < m; ++k) {
          B[r * m + k] -= f * B[c * m + k];
          Binv[r * m + k] -= f * Binv[c * m + k];
        }
      }
    }
    // x_B = -Binv * (sum over nonbasic j of A_j x_j)
    std::vector<double> rhs (m, 0.0);
    for (uint32_t j = 0; j < n; ++j) {
      if (rowOf[j] >= 0 || x[j] == 0.0) continue;
      forEachEntry (j, [&] (uint32_t i, double v) { rhs[i] += v * x[j]; });
    }
    for (size_t k = 0; k < m; ++k) {
      double s = 0.0;
      for (size_t i = 0; i < m; ++i) s -= Binv[k * m + i] * rhs[i];
      x[basis[k]] = s;
    }
    sinceRefactor = 0;
  }

  /** Duals of the current cost: y = c_B Binv. */
  void computeDuals (const std::vector<double> &c)
  {
    std::fill (y.begin (), y.end (), 0.0);
    for (size_t k = 0; k < m; ++k) {
      double ck = c[basis[k]];
      if (ck == 0.0) continue;
      const double *row = &Binv[k * m];
      for (size_t i = 0; i < m; ++i) y[i] += ck * row[i];
    }
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

  /** alpha = Binv * A_q. */
  void enteringColumn (uint32_t q)
  {
    std::fill (alpha.begin (), alpha.end (), 0.0);
    forEachEntry (q, [&] (uint32_t i, double v) {
      for (size_t k = 0; k < m; ++k) alpha[k] += Binv[k * m + i] * v;
    });
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
    // product form: row `leave` of Binv scaled, the others cleared of the entering column
    double piv = alpha[leave];
    double *prow = &Binv[leave * m];
    for (size_t i = 0; i < m; ++i) prow[i] /= piv;
    for (size_t k = 0; k < m; ++k) {
      if (k == leave || alpha[k] == 0.0) continue;
      double f = alpha[k];
      double *row = &Binv[k * m];
      for (size_t i = 0; i < m; ++i) row[i] -= f * prow[i];
    }
    rowOf[b] = -1;
    rowOf[q] = static_cast<int32_t> (leave);
    basis[leave] = q;
    ++pivots;
    ++sinceRefactor;
    if (sinceRefactor >= REFACTOR) refactor ();
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
    LpResult res;
    m = p.rows.size ();
    nStruct = p.columns;
    n = nStruct + 2 * m;
    if (m > limits.maxRows) {
      res.status = LpStatus::TooLarge;
      return res;
    }
    // columns from the rows
    cols.assign (nStruct, {});
    for (size_t i = 0; i < m; ++i)
      for (const auto &e : p.rows[i].coeffs) cols[e.first].emplace_back (static_cast<uint32_t> (i), static_cast<double> (e.second));
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
      lower[nStruct + i] = p.rows[i].lo;
      upper[nStruct + i] = p.rows[i].hi;
    }
    // row activities at the starting point, slacks clamped into their range, artificials take the rest
    std::vector<double> v (m, 0.0);
    for (size_t j = 0; j < nStruct; ++j)
      if (x[j] != 0.0)
        for (const auto &e : cols[j]) v[e.first] += e.second * x[j];
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
    Binv.assign (m * m, 0.0);
    for (size_t i = 0; i < m; ++i) Binv[i * m + i] = basis[i] < nStruct + m ? -1.0 : 1.0;
    y.assign (m, 0.0);
    alpha.assign (m, 0.0);
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
