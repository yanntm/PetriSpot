/*
 * CtlFormula.h
 *
 * The CTL formula tree: state predicates (expr::Expression) and the deadlock
 * atom at the leaves, booleans, and the path operators EX AX EF AF EG AG
 * EU AU EW AW (W is the weak until). Plain data; evaluation of the temporal
 * operators belongs to the checker (ctl/). See algorithm.md.
 */
#ifndef PETRI_EXPR_CTLFORMULA_H_
#define PETRI_EXPR_CTLFORMULA_H_

#include <algorithm>
#include <ostream>
#include <string>
#include <vector>

#include "expr/Expression.h"

namespace petri::expr
{

enum class CtlOp
{
  Pred,       // a state predicate, in pred
  Deadlock,   // the marking has no enabled transition
  NoDeadlock, // its negation
  Not, And, Or,
  EX, AX, EF, AF, EG, AG, // one child
  EU, AU, EW, AW          // two children: left, right
};

inline const char* to_string (CtlOp op)
{
  switch (op) {
  case CtlOp::Pred: return "pred";
  case CtlOp::Deadlock: return "deadlock";
  case CtlOp::NoDeadlock: return "(not deadlock)";
  case CtlOp::Not: return "not";
  case CtlOp::And: return "and";
  case CtlOp::Or: return "or";
  case CtlOp::EX: return "EX";
  case CtlOp::AX: return "AX";
  case CtlOp::EF: return "EF";
  case CtlOp::AF: return "AF";
  case CtlOp::EG: return "EG";
  case CtlOp::AG: return "AG";
  case CtlOp::EU: return "EU";
  case CtlOp::AU: return "AU";
  case CtlOp::EW: return "EW";
  case CtlOp::AW: return "AW";
  }
  return "?";
}

inline bool isTemporal (CtlOp op)
{
  return op >= CtlOp::EX;
}
inline bool isUnaryTemporal (CtlOp op)
{
  return op >= CtlOp::EX && op <= CtlOp::AG;
}
inline bool isBinaryTemporal (CtlOp op)
{
  return op >= CtlOp::EU;
}
/** Existential path quantifier. */
inline bool isExistential (CtlOp op)
{
  return op == CtlOp::EX || op == CtlOp::EF || op == CtlOp::EG || op == CtlOp::EU || op == CtlOp::EW;
}

struct CtlFormula
{
  CtlOp op = CtlOp::Pred;
  Expression pred; // Pred only
  std::vector<CtlFormula> kids;

  static CtlFormula predicate (Expression e)
  {
    CtlFormula f;
    f.op = CtlOp::Pred;
    f.pred = std::move (e);
    return f;
  }
  static CtlFormula constant (bool b)
  {
    return predicate (Expression::constant (b));
  }
  static CtlFormula leaf (CtlOp op)
  {
    CtlFormula f;
    f.op = op;
    return f;
  }
  static CtlFormula unary (CtlOp op, CtlFormula child)
  {
    CtlFormula f;
    f.op = op;
    f.kids.push_back (std::move (child));
    return f;
  }
  static CtlFormula binary (CtlOp op, CtlFormula left, CtlFormula right)
  {
    CtlFormula f;
    f.op = op;
    f.kids.push_back (std::move (left));
    f.kids.push_back (std::move (right));
    return f;
  }
  static CtlFormula nary (CtlOp op, std::vector<CtlFormula> children)
  {
    CtlFormula f;
    f.op = op;
    f.kids = std::move (children);
    return f;
  }

  /** A formula without temporal operator: a fact about one marking. */
  bool isState () const
  {
    if (isTemporal (op)) return false;
    for (const auto &k : kids) if (!k.isState ()) return false;
    return true;
  }
  bool isConstant () const
  {
    return op == CtlOp::Pred && pred.isConstant ();
  }
  bool isTrue () const
  {
    return op == CtlOp::Pred && pred.kind == Expression::Kind::True;
  }
  bool isFalse () const
  {
    return op == CtlOp::Pred && pred.kind == Expression::Kind::False;
  }
  const CtlFormula& left () const
  {
    return kids[0];
  }
  const CtlFormula& right () const
  {
    return kids.back ();
  }

  size_t size () const
  {
    size_t n = 1;
    for (const auto &k : kids) n += k.size ();
    return n;
  }
  /** Path quantifiers on the deepest branch. */
  size_t depth () const
  {
    size_t d = 0;
    for (const auto &k : kids) d = std::max (d, k.depth ());
    return d + (isTemporal (op) ? 1 : 0);
  }

  /**
   * Evaluate a state formula on a marking; deadlock is the flag given.
   * Undefined on a temporal formula.
   */
  template<class Marking>
    bool evalState (const Marking &m, bool deadlock) const
    {
      switch (op) {
      case CtlOp::Pred: return pred.eval (m);
      case CtlOp::Deadlock: return deadlock;
      case CtlOp::NoDeadlock: return !deadlock;
      case CtlOp::Not: return !kids[0].evalState (m, deadlock);
      case CtlOp::And:
        for (const auto &k : kids) if (!k.evalState (m, deadlock)) return false;
        return true;
      case CtlOp::Or:
        for (const auto &k : kids) if (k.evalState (m, deadlock)) return true;
        return false;
      default: return false;
      }
    }

  bool operator== (const CtlFormula &o) const
  {
    if (op != o.op) return false;
    if (op == CtlOp::Pred) return pred == o.pred;
    return kids == o.kids;
  }

  /** Infix form: E(a U b), AG(p >= 1), !(...), (a && b). */
  void print (std::ostream &os, const std::vector<std::string> *pnames = nullptr) const
  {
    switch (op) {
    case CtlOp::Pred: pred.print (os, pnames); break;
    case CtlOp::Deadlock: os << "deadlock"; break;
    case CtlOp::NoDeadlock: os << "!deadlock"; break;
    case CtlOp::Not:
      os << "!(";
      kids[0].print (os, pnames);
      os << ")";
      break;
    case CtlOp::And:
    case CtlOp::Or: {
      const char *sep = op == CtlOp::And ? " && " : " || ";
      os << "(";
      for (size_t i = 0; i < kids.size (); ++i) {
        if (i > 0) os << sep;
        kids[i].print (os, pnames);
      }
      os << ")";
      break;
    }
    case CtlOp::EX: case CtlOp::AX: case CtlOp::EF: case CtlOp::AF: case CtlOp::EG: case CtlOp::AG:
      os << to_string (op) << "(";
      kids[0].print (os, pnames);
      os << ")";
      break;
    case CtlOp::EU: case CtlOp::AU: case CtlOp::EW: case CtlOp::AW:
      os << (isExistential (op) ? "E(" : "A(");
      kids[0].print (os, pnames);
      os << (op == CtlOp::EU || op == CtlOp::AU ? " U " : " W ");
      kids[1].print (os, pnames);
      os << ")";
      break;
    }
  }
};

} // namespace petri::expr

#endif /* PETRI_EXPR_CTLFORMULA_H_ */
