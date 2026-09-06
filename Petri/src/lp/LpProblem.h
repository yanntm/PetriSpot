/*
 * LpProblem.h
 *
 * The statement of a linear program: columns with bounds, ranged rows of
 * integer coefficients (lo <= a.x <= hi, either side possibly infinite), an
 * optional linear objective to minimise. Plain data built by StateEquation
 * and read by Simplex; the integer coefficients are kept as such so an exact
 * checker can read the same problem the numerical solver saw. See
 * algorithm.md.
 */
#ifndef PETRI_LP_LPPROBLEM_H_
#define PETRI_LP_LPPROBLEM_H_

#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

namespace petri::lp
{

constexpr double INF = std::numeric_limits<double>::infinity ();

/** lo <= sum coeffs[i].second * x[coeffs[i].first] <= hi; coefficients sorted by column, non-zero. */
struct Row
{
  std::vector<std::pair<uint32_t, long long>> coeffs;
  double lo = -INF;
  double hi = INF;
};

struct LpProblem
{
  size_t columns = 0;
  std::vector<double> lower;    // per column, -INF allowed
  std::vector<double> upper;    // per column, INF allowed
  std::vector<Row> rows;
  std::vector<double> objective; // per column, to minimise; empty for a feasibility problem

  explicit LpProblem (size_t cols = 0)
      : columns (cols), lower (cols, 0.0), upper (cols, INF)
  {
  }

  bool hasObjective () const
  {
    return !objective.empty ();
  }

  /** Minimise sum of c[j] * x[j]. */
  void setObjective (std::vector<double> c)
  {
    objective = std::move (c);
  }

  /** Maximise the given form: stored negated, the solver minimises. */
  void setMaximise (const std::vector<double> &c)
  {
    objective.assign (c.size (), 0.0);
    for (size_t j = 0; j < c.size (); ++j) objective[j] = -c[j];
  }

  size_t addRow (Row r)
  {
    rows.push_back (std::move (r));
    return rows.size () - 1;
  }

  size_t nonZeros () const
  {
    size_t n = 0;
    for (const Row &r : rows) n += r.coeffs.size ();
    return n;
  }
};

} // namespace petri::lp

#endif /* PETRI_LP_LPPROBLEM_H_ */
