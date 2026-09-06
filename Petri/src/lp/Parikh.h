/*
 * Parikh.h
 *
 * A rational solution of the state equation into a Parikh hint: firing
 * counts rounded to the nearest integer, transitions with a positive count
 * only. The repair of a rounded vector that no longer keeps the marking
 * non-negative is planned (algorithm.md section 3); until then the walker
 * treats the hint as what it is, a suggestion.
 */
#ifndef PETRI_LP_PARIKH_H_
#define PETRI_LP_PARIKH_H_

#include <cmath>
#include <cstddef>
#include <vector>

#include "expr/Hint.h"

namespace petri::lp
{

inline petri::expr::ParikhHint toParikh (const std::vector<double> &x)
{
  petri::expr::ParikhHint h;
  for (size_t t = 0; t < x.size (); ++t) {
    long long k = std::llround (x[t]);
    if (k <= 0 && x[t] > 1e-6) k = 1; // a fractional firing is a firing
    if (k > 0) h.counts.emplace_back (t, k);
  }
  return h;
}

} // namespace petri::lp

#endif /* PETRI_LP_PARIKH_H_ */
