/*
 * Expression.h
 *
 * Plain-data boolean tree over linear atoms on place markings.
 * See algorithm.md for the shape and the normal form.
 */
#ifndef PETRI_EXPR_EXPRESSION_H_
#define PETRI_EXPR_EXPRESSION_H_

#include <algorithm>
#include <cstddef>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

namespace petri::expr
{

enum class Cmp
{
  EQ, NEQ, LE, GE, LT, GT
};

inline Cmp negate (Cmp c)
{
  switch (c) {
  case Cmp::EQ: return Cmp::NEQ;
  case Cmp::NEQ: return Cmp::EQ;
  case Cmp::LE: return Cmp::GT;
  case Cmp::GE: return Cmp::LT;
  case Cmp::LT: return Cmp::GE;
  case Cmp::GT: return Cmp::LE;
  }
  return c;
}

inline const char* to_string (Cmp c)
{
  switch (c) {
  case Cmp::EQ: return "==";
  case Cmp::NEQ: return "!=";
  case Cmp::LE: return "<=";
  case Cmp::GE: return ">=";
  case Cmp::LT: return "<";
  case Cmp::GT: return ">";
  }
  return "?";
}

inline bool compare (long long lhs, Cmp c, long long rhs)
{
  switch (c) {
  case Cmp::EQ: return lhs == rhs;
  case Cmp::NEQ: return lhs != rhs;
  case Cmp::LE: return lhs <= rhs;
  case Cmp::GE: return lhs >= rhs;
  case Cmp::LT: return lhs < rhs;
  case Cmp::GT: return lhs > rhs;
  }
  return false;
}

/**
 * sum_i coeff_i * m[place_i]  cmp  constant.
 * Terms are kept sorted by place with non-zero coefficients once normalised.
 */
struct LinearAtom
{
  std::vector<std::pair<size_t, long long>> terms;
  long long constant = 0;
  Cmp op = Cmp::LE;

  /** Accumulate coeff on place; call normalize() once all terms are added. */
  void addTerm (size_t place, long long coeff)
  {
    terms.emplace_back (place, coeff);
  }

  /** Sort by place, merge duplicates, drop zero coefficients. */
  void normalize ()
  {
    std::sort (terms.begin (), terms.end ());
    std::vector<std::pair<size_t, long long>> merged;
    merged.reserve (terms.size ());
    for (const auto &t : terms) {
      if (!merged.empty () && merged.back ().first == t.first) {
        merged.back ().second += t.second;
      } else {
        merged.push_back (t);
      }
    }
    merged.erase (
        std::remove_if (merged.begin (), merged.end (),
                        [] (const auto &t) { return t.second == 0; }),
        merged.end ());
    terms = std::move (merged);
  }

  template<class Marking>
  long long value (const Marking &m) const
  {
    long long v = 0;
    for (const auto &t : terms) {
      v += t.second * static_cast<long long> (m.get (t.first));
    }
    return v;
  }

  template<class Marking>
  bool holds (const Marking &m) const
  {
    return compare (value (m), op, constant);
  }

  bool operator== (const LinearAtom &o) const
  {
    return op == o.op && constant == o.constant && terms == o.terms;
  }

  void print (std::ostream &os, const std::vector<std::string> *pnames) const
  {
    if (terms.empty ()) os << "0";
    for (size_t i = 0; i < terms.size (); ++i) {
      long long c = terms[i].second;
      if (i > 0) os << (c < 0 ? " - " : " + ");
      else if (c < 0) os << "-";
      long long a = c < 0 ? -c : c;
      if (a != 1) os << a << "*";
      if (pnames) os << (*pnames)[terms[i].first];
      else os << "p" << terms[i].first;
    }
    os << " " << to_string (op) << " " << constant;
  }
};

struct Expression
{
  enum class Kind
  {
    True, False, Not, And, Or, Atom
  };

  Kind kind = Kind::True;
  std::vector<Expression> children;
  LinearAtom atom;

  static Expression constant (bool b)
  {
    Expression e;
    e.kind = b ? Kind::True : Kind::False;
    return e;
  }
  static Expression makeNot (Expression child)
  {
    Expression e;
    e.kind = Kind::Not;
    e.children.push_back (std::move (child));
    return e;
  }
  static Expression makeAnd (std::vector<Expression> kids)
  {
    Expression e;
    e.kind = Kind::And;
    e.children = std::move (kids);
    return e;
  }
  static Expression makeOr (std::vector<Expression> kids)
  {
    Expression e;
    e.kind = Kind::Or;
    e.children = std::move (kids);
    return e;
  }
  static Expression makeAtom (LinearAtom a)
  {
    Expression e;
    e.kind = Kind::Atom;
    e.atom = std::move (a);
    return e;
  }

  bool isConstant () const
  {
    return kind == Kind::True || kind == Kind::False;
  }

  template<class Marking>
  bool eval (const Marking &m) const
  {
    switch (kind) {
    case Kind::True: return true;
    case Kind::False: return false;
    case Kind::Not: return !children[0].eval (m);
    case Kind::And:
      for (const auto &c : children) if (!c.eval (m)) return false;
      return true;
    case Kind::Or:
      for (const auto &c : children) if (c.eval (m)) return true;
      return false;
    case Kind::Atom: return atom.holds (m);
    }
    return false;
  }

  /** Number of nodes in the tree. */
  size_t size () const
  {
    size_t n = 1;
    for (const auto &c : children) n += c.size ();
    return n;
  }

  bool operator== (const Expression &o) const
  {
    if (kind != o.kind) return false;
    if (kind == Kind::Atom) return atom == o.atom;
    return children == o.children;
  }

  void print (std::ostream &os,
              const std::vector<std::string> *pnames = nullptr) const
  {
    switch (kind) {
    case Kind::True: os << "true"; break;
    case Kind::False: os << "false"; break;
    case Kind::Not:
      os << "!(";
      children[0].print (os, pnames);
      os << ")";
      break;
    case Kind::And:
    case Kind::Or: {
      const char *sep = kind == Kind::And ? " && " : " || ";
      os << "(";
      for (size_t i = 0; i < children.size (); ++i) {
        if (i > 0) os << sep;
        children[i].print (os, pnames);
      }
      os << ")";
      break;
    }
    case Kind::Atom: atom.print (os, pnames); break;
    }
  }
};

} // namespace petri::expr

#endif /* PETRI_EXPR_EXPRESSION_H_ */
