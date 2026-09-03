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
  Unsupported   // parsed, outside the supported fragment; see comment
};

inline const char* to_string (PropertyKind k)
{
  switch (k) {
  case PropertyKind::Reachability: return "EF";
  case PropertyKind::Invariant: return "AG";
  case PropertyKind::Deadlock: return "EF deadlock";
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

  /** The state predicate the exploration has to reach (normal form). */
  Expression goal () const
  {
    if (kind == PropertyKind::Invariant) {
      return simplify (Expression::makeNot (body));
    }
    return simplify (body);
  }

  /** MCC verdict when a state satisfying goal() is reached. */
  const char* verdictIfReached () const
  {
    return kind == PropertyKind::Invariant ? "FALSE" : "TRUE";
  }

  void print (std::ostream &os,
              const std::vector<std::string> *pnames = nullptr) const
  {
    os << name << " [" << to_string (kind) << "]";
    if (kind == PropertyKind::Unsupported) {
      os << " : " << comment;
    } else if (kind != PropertyKind::Deadlock) {
      os << " : ";
      body.print (os, pnames);
    }
  }
};

} // namespace petri::expr

#endif /* PETRI_EXPR_PROPERTY_H_ */
