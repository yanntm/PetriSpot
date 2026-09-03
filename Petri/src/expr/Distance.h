/*
 * Distance.h
 *
 * Estimated distance from a marking to a normal-form Expression, in the
 * style of TAPAAL (Jensen, Nielsen, Oestergaard, Srba, ToPNoC 2016):
 * an atom's distance is how far its linear value is from satisfying the
 * comparison; And sums, Or takes the minimum. Zero iff the expression holds.
 */
#ifndef PETRI_EXPR_DISTANCE_H_
#define PETRI_EXPR_DISTANCE_H_

#include <cstdint>
#include <limits>

#include "expr/Expression.h"

namespace petri::expr
{

/** Distance of an unsatisfiable node; large but safe to add a few times. */
constexpr uint64_t INFINITE_DISTANCE = std::numeric_limits<uint64_t>::max () / 1024;

inline uint64_t atomDistance (long long value, Cmp op, long long k)
{
  switch (op) {
  case Cmp::LE:
  case Cmp::LT:
    return value <= k ? 0 : static_cast<uint64_t> (value - k);
  case Cmp::GE:
  case Cmp::GT:
    return value >= k ? 0 : static_cast<uint64_t> (k - value);
  case Cmp::EQ:
    return value == k ? 0 : static_cast<uint64_t> (value > k ? value - k : k - value);
  case Cmp::NEQ:
    return value != k ? 0 : 1;
  }
  return 0;
}

template<class Marking>
  uint64_t distance (const Expression &e, const Marking &m)
  {
    switch (e.kind) {
    case Expression::Kind::True: return 0;
    case Expression::Kind::False: return INFINITE_DISTANCE;
    case Expression::Kind::Not: {
      // not in normal form; fall back to exact evaluation
      return e.eval (m) ? 0 : 1;
    }
    case Expression::Kind::And: {
      uint64_t sum = 0;
      for (const auto &c : e.children) {
        uint64_t d = distance (c, m);
        if (d >= INFINITE_DISTANCE) return INFINITE_DISTANCE;
        sum += d;
      }
      return sum;
    }
    case Expression::Kind::Or: {
      uint64_t best = INFINITE_DISTANCE;
      for (const auto &c : e.children) {
        uint64_t d = distance (c, m);
        if (d < best) best = d;
        if (best == 0) break;
      }
      return best;
    }
    case Expression::Kind::Atom:
      return atomDistance (e.atom.value (m), e.atom.op, e.atom.constant);
    }
    return INFINITE_DISTANCE;
  }

} // namespace petri::expr

#endif /* PETRI_EXPR_DISTANCE_H_ */
