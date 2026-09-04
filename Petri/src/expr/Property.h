/*
 * Property.h
 *
 * A named property over a net: kind, body, and the exploration goal derived
 * from them. See algorithm.md.
 */
#ifndef PETRI_EXPR_PROPERTY_H_
#define PETRI_EXPR_PROPERTY_H_

#include <string>
#include <vector>

#include "expr/Expression.h"
#include "expr/Simplify.h"

namespace petri::expr
{

enum class PropertyKind
{
  Reachability, // EF body
  Invariant,    // AG body
  Deadlock,     // EF deadlock
  Bound,        // maximise body (a single linear atom holds the form)
  Unsupported   // parsed, outside the supported fragment; see comment
};

inline const char* to_string (PropertyKind k)
{
  switch (k) {
  case PropertyKind::Reachability: return "EF";
  case PropertyKind::Invariant: return "AG";
  case PropertyKind::Deadlock: return "EF deadlock";
  case PropertyKind::Bound: return "bound";
  case PropertyKind::Unsupported: return "unsupported";
  }
  return "?";
}

struct Property
{
  std::string name;
  PropertyKind kind = PropertyKind::Unsupported;
  Expression body;
  std::string comment;
  long long boundHint = -1; // Bound: a known upper bound of the form, -1 when unknown

  /** Bound: the linear form to maximise (the terms of body's atom). */
  const LinearAtom& boundForm () const
  {
    return body.atom;
  }

  /**
   * The state predicate the exploration has to reach (normal form). For a
   * bound with a hint it is "form >= hint"; a bound without hint has no
   * reachable goal (False).
   */
  Expression goal () const
  {
    if (kind == PropertyKind::Invariant) {
      return simplify (Expression::makeNot (body));
    }
    if (kind == PropertyKind::Bound) {
      if (boundHint < 0) return Expression::constant (false);
      LinearAtom a = body.atom;
      a.op = Cmp::GE;
      a.constant = boundHint;
      return simplify (Expression::makeAtom (std::move (a)));
    }
    return simplify (body);
  }

  /** MCC verdict when a state satisfying goal() is reached. */
  std::string verdictIfReached () const
  {
    if (kind == PropertyKind::Bound) return std::to_string (boundHint);
    return kind == PropertyKind::Invariant ? "FALSE" : "TRUE";
  }

  void print (std::ostream &os,
              const std::vector<std::string> *pnames = nullptr) const
  {
    os << name << " [" << to_string (kind) << "]";
    if (kind == PropertyKind::Unsupported) {
      os << " : " << comment;
    } else if (kind == PropertyKind::Bound) {
      os << " : max ";
      LinearAtom form = body.atom;
      form.op = Cmp::LE;
      form.constant = boundHint;
      form.print (os, pnames);
      os << (boundHint < 0 ? " (no hint)" : " (hint)");
    } else if (kind != PropertyKind::Deadlock) {
      os << " : ";
      body.print (os, pnames);
    }
  }
};

} // namespace petri::expr

#endif /* PETRI_EXPR_PROPERTY_H_ */
