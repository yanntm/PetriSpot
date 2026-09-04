/*
 * Simplify.h
 *
 * Normal form of an Expression: no Not nodes, integer comparisons restricted
 * to ==, !=, <=, >=, constants folded, And/Or flattened and deduplicated.
 * See algorithm.md.
 */
#ifndef PETRI_EXPR_SIMPLIFY_H_
#define PETRI_EXPR_SIMPLIFY_H_

#include "expr/Expression.h"

namespace petri::expr
{

namespace detail
{

/**
 * Normalise an atom: sorted merged terms, strict comparisons rewritten,
 * leading coefficient positive (the atom is negated as a whole when needed),
 * and folding when no place remains or when the non-negativity of markings
 * decides the comparison (all coefficients positive, constant out of range).
 */
inline Expression simplifyAtom (LinearAtom a)
{
  a.normalize ();
  if (a.op == Cmp::LT) {
    a.op = Cmp::LE;
    a.constant -= 1;
  } else if (a.op == Cmp::GT) {
    a.op = Cmp::GE;
    a.constant += 1;
  }
  if (a.terms.empty ()) {
    return Expression::constant (compare (0, a.op, a.constant));
  }
  if (a.terms[0].second < 0) {
    for (auto &t : a.terms) t.second = -t.second;
    a.constant = -a.constant;
    if (a.op == Cmp::LE) a.op = Cmp::GE;
    else if (a.op == Cmp::GE) a.op = Cmp::LE;
  }
  bool allPositive = true;
  for (const auto &t : a.terms) if (t.second < 0) allPositive = false;
  if (allPositive) {
    // the value is a non-negative integer
    switch (a.op) {
    case Cmp::LE: if (a.constant < 0) return Expression::constant (false); break;
    case Cmp::GE: if (a.constant <= 0) return Expression::constant (true); break;
    case Cmp::EQ: if (a.constant < 0) return Expression::constant (false); break;
    case Cmp::NEQ: if (a.constant < 0) return Expression::constant (true); break;
    default: break;
    }
  }
  return Expression::makeAtom (std::move (a));
}

Expression simplifyRec (const Expression &e, bool negated);

/** Build a flattened, deduplicated And (isAnd) or Or over simplified kids. */
inline Expression simplifyNary (const Expression &e, bool isAnd, bool negated)
{
  const Expression::Kind kind = isAnd ? Expression::Kind::And : Expression::Kind::Or;
  const Expression::Kind neutral = isAnd ? Expression::Kind::True : Expression::Kind::False;
  const Expression::Kind absorbing = isAnd ? Expression::Kind::False : Expression::Kind::True;
  std::vector<Expression> kids;
  for (const auto &c : e.children) {
    Expression s = simplifyRec (c, negated);
    if (s.kind == absorbing) return Expression::constant (!isAnd);
    if (s.kind == neutral) continue;
    if (s.kind == kind) {
      for (auto &g : s.children) kids.push_back (std::move (g));
    } else {
      kids.push_back (std::move (s));
    }
  }
  std::vector<Expression> unique;
  for (auto &k : kids) {
    bool seen = false;
    for (const auto &u : unique) {
      if (u == k) { seen = true; break; }
    }
    if (!seen) unique.push_back (std::move (k));
  }
  if (unique.empty ()) return Expression::constant (isAnd);
  if (unique.size () == 1) return std::move (unique[0]);
  return isAnd ? Expression::makeAnd (std::move (unique)) : Expression::makeOr (std::move (unique));
}

inline Expression simplifyRec (const Expression &e, bool negated)
{
  switch (e.kind) {
  case Expression::Kind::True: return Expression::constant (!negated);
  case Expression::Kind::False: return Expression::constant (negated);
  case Expression::Kind::Not: return simplifyRec (e.children[0], !negated);
  case Expression::Kind::And:
    // !(a && b) == (!a || !b): the negated conjunction is a disjunction
    return simplifyNary (e, !negated, negated);
  case Expression::Kind::Or:
    return simplifyNary (e, negated, negated);
  case Expression::Kind::Atom: {
    LinearAtom a = e.atom;
    if (negated) a.op = negate (a.op);
    return simplifyAtom (std::move (a));
  }
  }
  return Expression::constant (false);
}

} // namespace detail

/** Equivalent expression in the normal form of algorithm.md. */
inline Expression simplify (const Expression &e)
{
  return detail::simplifyRec (e, false);
}

} // namespace petri::expr

#endif /* PETRI_EXPR_SIMPLIFY_H_ */
