/*
 * Basis.h
 *
 * The inverse of the simplex basis in product form: a signed diagonal (the
 * slack and artificial columns of the start) followed by one eta matrix per
 * pivot, each the identity except for one column, kept as the sparse pivot
 * column it came from. The inverse of a sparse basis is dense and is never
 * formed; it is applied. ftran solves B v' = v (an entering column, the
 * basic values), btran solves y' B = y (the duals). Work and storage are the
 * total size of the eta columns; a refactorisation rebuilds the file from
 * the basis columns when it has grown long. See algorithm.md section 3.
 */
#ifndef PETRI_LP_BASIS_H_
#define PETRI_LP_BASIS_H_

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace petri::lp
{

class Basis
{
  struct Eta
  {
    uint32_t row;    // the pivot row r
    double pivot;    // alpha_r
    std::vector<std::pair<uint32_t, double>> col; // alpha_k for k != r, non-zero
  };
  static constexpr double DROP = 1e-14; // an eta entry below this is not stored

  size_t m = 0;
  std::vector<double> diag;  // the start: B0 = diag, so B0^-1 = diag (entries +1 or -1)
  std::vector<Eta> etas;
  size_t stored = 0;         // entries in the eta file

public:
  /** Start over from a signed diagonal basis. */
  void reset (size_t rows, const std::vector<double> &diagonal)
  {
    m = rows;
    diag = diagonal;
    etas.clear ();
    stored = 0;
  }

  size_t etaCount () const
  {
    return etas.size ();
  }
  size_t nonZeros () const
  {
    return stored;
  }

  /** v <- B^-1 v, in place; v is dense of size m. */
  void ftran (std::vector<double> &v) const
  {
    for (size_t i = 0; i < m; ++i) v[i] *= diag[i];
    for (const Eta &e : etas) {
      double vr = v[e.row] / e.pivot;
      v[e.row] = vr;
      if (vr == 0.0) continue;
      for (const auto &a : e.col) v[a.first] -= a.second * vr;
    }
  }

  /** y <- y B^-1 (y a row vector), in place; y is dense of size m. */
  void btran (std::vector<double> &y) const
  {
    for (size_t k = etas.size (); k-- > 0;) {
      const Eta &e = etas[k];
      double s = y[e.row];
      for (const auto &a : e.col) s -= y[a.first] * a.second;
      y[e.row] = s / e.pivot;
    }
    for (size_t i = 0; i < m; ++i) y[i] *= diag[i];
  }

  /** Record the pivot at `row` of the column alpha = B^-1 A_q, so that B^-1 now has A_q in that row's basic slot. */
  void pivot (size_t row, const std::vector<double> &alpha)
  {
    Eta e;
    e.row = static_cast<uint32_t> (row);
    e.pivot = alpha[row];
    for (size_t k = 0; k < m; ++k) {
      if (k == row || std::fabs (alpha[k]) < DROP) continue;
      e.col.emplace_back (static_cast<uint32_t> (k), alpha[k]);
    }
    stored += e.col.size () + 1;
    etas.push_back (std::move (e));
  }
};

} // namespace petri::lp

#endif /* PETRI_LP_BASIS_H_ */
