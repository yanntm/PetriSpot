/*
 * CtlSimplify.h
 *
 * Negation normal form and simplification of CTL formulas: negation pushed
 * to the leaves through the duals of the operators, constants folded,
 * booleans over state leaves merged into one predicate, and the collapsing
 * rewrite rules of nested operators (EF EF a = EF a, EF AF a = EF a,
 * A[a U EF b] = EF b, ...). Also the two state-predicate approximations of
 * a formula used to steer a walk: now(f), implied by f, and suf(f),
 * implying f. See algorithm.md.
 */
#ifndef PETRI_EXPR_CTLSIMPLIFY_H_
#define PETRI_EXPR_CTLSIMPLIFY_H_

#include <utility>
#include <vector>

#include "expr/CtlFormula.h"
#include "expr/Simplify.h"

namespace petri::expr
{

CtlFormula ctlNnf (const CtlFormula &f);

/** The negation of f in negation normal form. */
inline CtlFormula ctlNegate (const CtlFormula &f)
{
  using Op = CtlOp;
  switch (f.op) {
  case Op::Pred: return CtlFormula::predicate (simplify (Expression::makeNot (f.pred)));
  case Op::Deadlock: return CtlFormula::leaf (Op::NoDeadlock);
  case Op::NoDeadlock: return CtlFormula::leaf (Op::Deadlock);
  case Op::Not: return ctlNnf (f.kids[0]);
  case Op::And:
  case Op::Or: {
    std::vector<CtlFormula> ks;
    for (const auto &k : f.kids) ks.push_back (ctlNegate (k));
    return CtlFormula::nary (f.op == Op::And ? Op::Or : Op::And, std::move (ks));
  }
  case Op::EX: return CtlFormula::unary (Op::AX, ctlNegate (f.kids[0]));
  case Op::AX: return CtlFormula::unary (Op::EX, ctlNegate (f.kids[0]));
  case Op::EF: return CtlFormula::unary (Op::AG, ctlNegate (f.kids[0]));
  case Op::AG: return CtlFormula::unary (Op::EF, ctlNegate (f.kids[0]));
  case Op::EG: return CtlFormula::unary (Op::AF, ctlNegate (f.kids[0]));
  case Op::AF: return CtlFormula::unary (Op::EG, ctlNegate (f.kids[0]));
  case Op::EU: case Op::AU: case Op::EW: case Op::AW: {
    // not E[a U b] = A[not b W (not a and not b)], not E[a W b] = A[not b U (not a and not b)]
    CtlFormula nb = ctlNegate (f.kids[1]);
    CtlFormula both = CtlFormula::nary (Op::And, { ctlNegate (f.kids[0]), nb });
    Op d = f.op == Op::EU ? Op::AW : f.op == Op::AU ? Op::EW : f.op == Op::EW ? Op::AU : Op::EU;
    return CtlFormula::binary (d, std::move (nb), std::move (both));
  }
  }
  return f;
}

/** Negation normal form: Not only remains inside state predicates, folded by simplify. */
inline CtlFormula ctlNnf (const CtlFormula &f)
{
  if (f.op == CtlOp::Not) return ctlNegate (f.kids[0]);
  if (f.op == CtlOp::Pred) return CtlFormula::predicate (simplify (f.pred));
  CtlFormula r;
  r.op = f.op;
  for (const auto &k : f.kids) r.kids.push_back (ctlNnf (k));
  return r;
}

CtlFormula ctlSimplify (const CtlFormula &f);

namespace detail
{

/** And/Or: flatten, merge the predicate children into one, fold constants, single child. */
inline CtlFormula simplifyBool (CtlOp op, std::vector<CtlFormula> kids)
{
  bool isAnd = op == CtlOp::And;
  std::vector<CtlFormula> flat;
  for (auto &k : kids) {
    if (k.op == op) for (auto &kk : k.kids) flat.push_back (std::move (kk));
    else flat.push_back (std::move (k));
  }
  std::vector<Expression> preds;
  std::vector<CtlFormula> rest;
  for (auto &k : flat) {
    if (k.op == CtlOp::Pred) preds.push_back (std::move (k.pred));
    else {
      bool dup = false;
      for (const auto &r : rest) if (r == k) { dup = true; break; }
      if (!dup) rest.push_back (std::move (k));
    }
  }
  if (!preds.empty ()) {
    Expression merged = simplify (isAnd ? Expression::makeAnd (std::move (preds)) : Expression::makeOr (std::move (preds)));
    bool absorbing = isAnd ? merged.kind == Expression::Kind::False : merged.kind == Expression::Kind::True;
    bool neutral = isAnd ? merged.kind == Expression::Kind::True : merged.kind == Expression::Kind::False;
    if (absorbing) return CtlFormula::constant (!isAnd);
    if (!neutral) rest.insert (rest.begin (), CtlFormula::predicate (std::move (merged)));
  }
  // deadlock beside its negation
  bool dl = false, ndl = false;
  for (const auto &r : rest) {
    if (r.op == CtlOp::Deadlock) dl = true;
    if (r.op == CtlOp::NoDeadlock) ndl = true;
  }
  if (dl && ndl) return CtlFormula::constant (!isAnd);
  if (rest.empty ()) return CtlFormula::constant (isAnd);
  if (rest.size () == 1) return std::move (rest[0]);
  return CtlFormula::nary (op, std::move (rest));
}

/** Or(ks) with EF children pulled out: the EF parts and the others, separately. */
inline void splitEF (const CtlFormula &g, std::vector<CtlFormula> &efs, std::vector<CtlFormula> &others)
{
  if (g.op == CtlOp::Or) {
    for (const auto &k : g.kids) (k.op == CtlOp::EF ? efs : others).push_back (k);
  } else if (g.op == CtlOp::EF) {
    efs.push_back (g);
  } else {
    others.push_back (g);
  }
}

inline CtlFormula simplifyUnary (CtlOp op, CtlFormula g)
{
  using Op = CtlOp;
  if (g.isConstant ()) {
    // F and G of a constant are that constant; X of a constant speaks of the successors' existence
    if (op == Op::EX) return g.isTrue () ? CtlFormula::leaf (Op::NoDeadlock) : g;
    if (op == Op::AX) return g.isFalse () ? CtlFormula::leaf (Op::Deadlock) : g;
    return g;
  }
  switch (op) {
  case Op::EF:
    if (g.op == Op::EF || g.op == Op::AF) return g.op == Op::EF ? g : CtlFormula::unary (Op::EF, g.kids[0]);
    if (g.op == Op::EU || g.op == Op::AU) return ctlSimplify (CtlFormula::unary (Op::EF, g.kids[1]));
    if (g.op == Op::Or) {
      std::vector<CtlFormula> ks;
      for (auto &k : g.kids) ks.push_back (ctlSimplify (CtlFormula::unary (Op::EF, std::move (k))));
      return simplifyBool (Op::Or, std::move (ks));
    }
    break;
  case Op::AF:
    if (g.op == Op::EF || g.op == Op::AF) return g;
    if (g.op == Op::AU) return ctlSimplify (CtlFormula::unary (Op::AF, g.kids[1]));
    if (g.op == Op::Or) {
      std::vector<CtlFormula> efs, others;
      splitEF (g, efs, others);
      if (!efs.empty ()) {
        if (!others.empty ()) efs.push_back (CtlFormula::unary (Op::AF, simplifyBool (Op::Or, std::move (others))));
        return simplifyBool (Op::Or, std::move (efs));
      }
    }
    break;
  case Op::EG:
    if (g.op == Op::EG || g.op == Op::AG) return g;
    break;
  case Op::AG:
    if (g.op == Op::AG) return g;
    if (g.op == Op::And) {
      std::vector<CtlFormula> ks;
      for (auto &k : g.kids) ks.push_back (ctlSimplify (CtlFormula::unary (Op::AG, std::move (k))));
      return simplifyBool (Op::And, std::move (ks));
    }
    break;
  default: break;
  }
  return CtlFormula::unary (op, std::move (g));
}

inline CtlFormula simplifyBinary (CtlOp op, CtlFormula a, CtlFormula b)
{
  using Op = CtlOp;
  bool until = op == Op::EU || op == Op::AU;
  bool exist = op == Op::EU || op == Op::EW;
  if (b.isTrue ()) return b;
  if (a.isFalse ()) return b; // must hold now
  if (until) {
    if (b.isFalse ()) return b;
    if (a.isTrue ()) return ctlSimplify (CtlFormula::unary (exist ? Op::EF : Op::AF, std::move (b)));
    if (a.op == Op::Deadlock) return b;      // stuck at s: b now or never
    if (a.op == Op::NoDeadlock) return ctlSimplify (CtlFormula::unary (exist ? Op::EF : Op::AF, std::move (b)));
    if (b.op == Op::NoDeadlock) return b;    // holds now, or s is a deadlock and nothing moves
    if (b.op == Op::EF) return b;            // X[a U EF c] = EF c
    if (!exist && b.op == Op::AF) return b;  // A[a U AF c] = AF c
    std::vector<CtlFormula> efs, others;
    splitEF (b, efs, others);
    if (!efs.empty () && b.op == Op::Or) {
      if (!others.empty ()) efs.push_back (CtlFormula::binary (op, a, simplifyBool (Op::Or, std::move (others))));
      return simplifyBool (Op::Or, std::move (efs));
    }
  } else {
    if (a.isTrue ()) return a;
    if (b.isFalse ()) return ctlSimplify (CtlFormula::unary (exist ? Op::EG : Op::AG, std::move (a)));
  }
  return CtlFormula::binary (op, std::move (a), std::move (b));
}

} // namespace detail

/** Simplify a formula in negation normal form (apply ctlNnf first). */
inline CtlFormula ctlSimplify (const CtlFormula &f)
{
  using Op = CtlOp;
  switch (f.op) {
  case Op::Pred: return CtlFormula::predicate (simplify (f.pred));
  case Op::Deadlock:
  case Op::NoDeadlock: return f;
  case Op::Not: return ctlSimplify (ctlNegate (f.kids[0]));
  case Op::And:
  case Op::Or: {
    std::vector<CtlFormula> ks;
    for (const auto &k : f.kids) ks.push_back (ctlSimplify (k));
    return detail::simplifyBool (f.op, std::move (ks));
  }
  case Op::EX: case Op::AX: case Op::EF: case Op::AF: case Op::EG: case Op::AG:
    return detail::simplifyUnary (f.op, ctlSimplify (f.kids[0]));
  case Op::EU: case Op::AU: case Op::EW: case Op::AW:
    return detail::simplifyBinary (f.op, ctlSimplify (f.kids[0]), ctlSimplify (f.kids[1]));
  }
  return f;
}

/** Normal form of any formula: negation pushed down, then simplified. */
inline CtlFormula ctlNormalize (const CtlFormula &f)
{
  return ctlSimplify (ctlNnf (f));
}

/** The negation of a normal-form formula, in normal form. */
inline CtlFormula ctlDual (const CtlFormula &f)
{
  return ctlSimplify (ctlNegate (f));
}

/**
 * A state predicate implied by f at any state where f holds (weak: true
 * for X and F nodes, whose truth says nothing about the current state).
 */
inline Expression ctlNow (const CtlFormula &f)
{
  using Op = CtlOp;
  switch (f.op) {
  case Op::Pred: return f.pred;
  case Op::And:
  case Op::Or: {
    std::vector<Expression> ks;
    for (const auto &k : f.kids) ks.push_back (ctlNow (k));
    return simplify (f.op == Op::And ? Expression::makeAnd (std::move (ks)) : Expression::makeOr (std::move (ks)));
  }
  case Op::EG: case Op::AG: return ctlNow (f.kids[0]);
  case Op::EU: case Op::AU: case Op::EW: case Op::AW:
    return simplify (Expression::makeOr ({ ctlNow (f.kids[0]), ctlNow (f.kids[1]) }));
  default: return Expression::constant (true);
  }
}

/**
 * A state predicate implying f wherever it holds (false when there is none:
 * X and G nodes, the deadlock leaves).
 */
inline Expression ctlSuf (const CtlFormula &f)
{
  using Op = CtlOp;
  switch (f.op) {
  case Op::Pred: return f.pred;
  case Op::And:
  case Op::Or: {
    std::vector<Expression> ks;
    for (const auto &k : f.kids) ks.push_back (ctlSuf (k));
    return simplify (f.op == Op::And ? Expression::makeAnd (std::move (ks)) : Expression::makeOr (std::move (ks)));
  }
  case Op::EF: case Op::AF: return ctlSuf (f.kids[0]);
  case Op::EU: case Op::AU: case Op::EW: case Op::AW: return ctlSuf (f.kids[1]);
  default: return Expression::constant (false);
  }
}

} // namespace petri::expr

#endif /* PETRI_EXPR_CTLSIMPLIFY_H_ */
