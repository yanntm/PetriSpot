/*
 * InitialState.h
 *
 * What a property says about the initial marking alone, before any
 * exploration: the three-valued truth of a CTL formula at the initial state,
 * the rewrite of an until whose left side already fails there, and the
 * requalification of a formula whose shape after that is the reachability
 * fragment. The rules follow Bonneland et al., "Simplification of CTL
 * Formulae for Efficient Model Checking of Petri Nets" (Petri Nets 2018),
 * section 3, as ITS-Tools applies them. See algorithm.md, "Initial state".
 */
#ifndef PETRI_EXPR_INITIALSTATE_H_
#define PETRI_EXPR_INITIALSTATE_H_

#include <vector>

#include "core/SparsePetriNet.h"
#include "expr/CtlSimplify.h"
#include "expr/Property.h"

namespace petri::expr
{

enum class Tri { False, Unknown, True };

inline Tri triOf (bool b) { return b ? Tri::True : Tri::False; }
inline Tri triNot (Tri v) { return v == Tri::True ? Tri::False : v == Tri::False ? Tri::True : Tri::Unknown; }

/** Whether no transition of the net is enabled at its initial marking. */
template<class T>
  bool initialDeadlock (const SparsePetriNet<T> &net)
  {
    const auto &marks = net.getMarks ();
    for (size_t t = 0; t < net.getTransitionCount (); ++t) {
      const auto &pre = net.getFlowPT ().getColumn (t);
      bool enabled = true;
      for (size_t i = 0; i < pre.size () && enabled; ++i) enabled = marks[pre.keyAt (i)] >= pre.valueAt (i);
      if (enabled) return false;
    }
    return true;
  }

/**
 * Truth of f at the initial marking m, when its shape decides it there:
 * a state formula is evaluated; G of a false operand is false, F of a true
 * operand is true; an until is true when its right side holds and false
 * when neither side holds; X says nothing. Unknown otherwise.
 */
template<class Marking>
  Tri initialValue (const CtlFormula &f, const Marking &m, bool deadlock)
  {
    using Op = CtlOp;
    switch (f.op) {
    case Op::Pred: return triOf (f.pred.eval (m));
    case Op::Deadlock: return triOf (deadlock);
    case Op::NoDeadlock: return triOf (!deadlock);
    case Op::Not: return triNot (initialValue (f.kids[0], m, deadlock));
    case Op::And:
    case Op::Or: {
      bool isAnd = f.op == Op::And;
      Tri absorbing = isAnd ? Tri::False : Tri::True, neutral = triNot (absorbing);
      bool unknown = false;
      for (const auto &k : f.kids) {
        Tri v = initialValue (k, m, deadlock);
        if (v == absorbing) return absorbing;
        if (v == Tri::Unknown) unknown = true;
      }
      return unknown ? Tri::Unknown : neutral;
    }
    case Op::EX:
    case Op::AX: return Tri::Unknown;
    case Op::EG:
    case Op::AG: return initialValue (f.kids[0], m, deadlock) == Tri::False ? Tri::False : Tri::Unknown;
    case Op::EF:
    case Op::AF: return initialValue (f.kids[0], m, deadlock) == Tri::True ? Tri::True : Tri::Unknown;
    case Op::EU: case Op::AU: case Op::EW: case Op::AW: {
      Tri right = initialValue (f.kids[1], m, deadlock);
      if (right == Tri::True) return Tri::True;
      if (right == Tri::False && initialValue (f.kids[0], m, deadlock) == Tri::False) return Tri::False;
      return Tri::Unknown;
    }
    }
    return Tri::Unknown;
  }

/**
 * f rewritten with what the initial marking decides, for a formula in
 * negation normal form. A decided formula becomes its constant; the children
 * of a boolean are rewritten in turn (they too speak of the initial state);
 * an until whose left side fails initially is its right side, which then
 * has to hold at the initial state. Operands under a path operator are left
 * alone: they speak of other states.
 */
template<class Marking>
  CtlFormula initialRewrite (const CtlFormula &f, const Marking &m, bool deadlock)
  {
    using Op = CtlOp;
    Tri v = initialValue (f, m, deadlock);
    if (v != Tri::Unknown) return CtlFormula::constant (v == Tri::True);
    switch (f.op) {
    case Op::Not: return ctlSimplify (ctlNegate (initialRewrite (f.kids[0], m, deadlock)));
    case Op::And:
    case Op::Or: {
      std::vector<CtlFormula> ks;
      for (const auto &k : f.kids) ks.push_back (initialRewrite (k, m, deadlock));
      return ctlSimplify (CtlFormula::nary (f.op, std::move (ks)));
    }
    case Op::EU: case Op::AU: case Op::EW: case Op::AW:
      if (initialValue (f.kids[0], m, deadlock) == Tri::False) return initialRewrite (f.kids[1], m, deadlock);
      return f;
    default: return f;
    }
  }

/**
 * Give a CTL property the kind its normal form deserves: EF of a state
 * predicate is a reachability property, AG of one an invariant, EF deadlock
 * the deadlock property. The formula is kept as CTL otherwise.
 */
inline void requalify (Property &p)
{
  if (p.kind != PropertyKind::CTL) return;
  const CtlFormula &f = p.ctl;
  if (f.op == CtlOp::EF && f.kids[0].op == CtlOp::Pred) {
    p.kind = PropertyKind::Reachability;
    p.body = f.kids[0].pred;
  } else if (f.op == CtlOp::AG && f.kids[0].op == CtlOp::Pred) {
    p.kind = PropertyKind::Invariant;
    p.body = f.kids[0].pred;
  } else if (f.op == CtlOp::EF && f.kids[0].op == CtlOp::Deadlock) {
    p.kind = PropertyKind::Deadlock;
  } else {
    return;
  }
  p.ctl = CtlFormula ();
}

/**
 * Simplify every property with the initial marking of the net: a
 * reachability or invariant body decided there becomes a constant, the
 * deadlock property is TRUE on an initially dead net, a CTL formula is
 * rewritten as above then requalified. Returns how many properties changed.
 */
template<class T>
  size_t simplifyInitial (const SparsePetriNet<T> &net, std::vector<Property> &properties)
  {
    const SparseArray<T> initial (net.getMarks ());
    bool deadlock = initialDeadlock (net);
    size_t changed = 0;
    for (auto &p : properties) {
      switch (p.kind) {
      case PropertyKind::Reachability:
        if (!p.body.isConstant () && p.body.eval (initial)) { p.body = Expression::constant (true); ++changed; }
        break;
      case PropertyKind::Invariant:
        if (!p.body.isConstant () && !p.body.eval (initial)) { p.body = Expression::constant (false); ++changed; }
        break;
      case PropertyKind::Deadlock:
        if (deadlock) { p.kind = PropertyKind::Reachability; p.body = Expression::constant (true); ++changed; }
        break;
      case PropertyKind::CTL: {
        CtlFormula f = ctlNormalize (p.ctl);
        f = initialRewrite (f, initial, deadlock);
        if (!(f == p.ctl)) ++changed;
        p.ctl = std::move (f);
        requalify (p);
        break;
      }
      default: break;
      }
    }
    return changed;
  }

} // namespace petri::expr

#endif /* PETRI_EXPR_INITIALSTATE_H_ */
