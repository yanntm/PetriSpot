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
#include <limits>
#include <utility>
#include <vector>

#include "core/MatrixCol.h"
#include "core/SparseArray.h"

namespace petri::lp
{

constexpr double INF = std::numeric_limits<double>::infinity ();

/** lo <= coeffs . x <= hi, coeffs a sparse vector over the columns. */
struct Row
{
  SparseArray<long long> coeffs;
  double lo = -INF;
  double hi = INF;
};

struct LpProblem
{
  size_t columns = 0;
  std::vector<double> lower;    // per column, -INF allowed
  std::vector<double> upper;    // per column, INF allowed
  MatrixCol<long long> rows;    // one sparse column of this matrix per row of the program, over the columns
  std::vector<double> lo, hi;   // per row
  std::vector<double> objective; // per column, to minimise; empty for a feasibility problem

  explicit LpProblem (size_t cols = 0)
      : columns (cols), lower (cols, 0.0), upper (cols, INF), rows (cols, 0)
  {
  }

  size_t rowCount () const
  {
    return rows.getColumnCount ();
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
    colsValid = false;
    rows.appendColumn (std::move (r.coeffs));
    lo.push_back (r.lo);
    hi.push_back (r.hi);
    return rowCount () - 1;
  }

  /** The program's matrix by column: the transpose of the rows, computed once and kept while no row is added. */
  const MatrixCol<long long>& byColumn () const
  {
    if (!colsValid) {
      cols = rows.transpose ();
      colsValid = true;
    }
    return cols;
  }

private:
  mutable MatrixCol<long long> cols;
  mutable bool colsValid = false;
};

} // namespace petri::lp

#endif /* PETRI_LP_LPPROBLEM_H_ */
