/*
 * SexprPrinter.h
 *
 * Print expressions and properties in the s-expression syntax read by
 * parse/sexpr/ (INTEROP.md section 4): places as p<i> / transitions as t<i>
 * when no name table is given, quoted names otherwise when needed.
 */
#ifndef PETRI_EXPR_SEXPRPRINTER_H_
#define PETRI_EXPR_SEXPRPRINTER_H_

#include <cctype>
#include <ostream>
#include <string>
#include <vector>

#include "expr/Expression.h"
#include "expr/Property.h"

namespace petri::expr
{

/** A name is printed bare when the reader would take it back as the same name. */
inline bool needsQuotes (const std::string &s)
{
  if (s.empty ()) return true;
  bool indexLike = (s[0] == 'p' || s[0] == 't') && s.size () > 1;
  for (size_t i = 0; i < s.size (); ++i) {
    unsigned char c = static_cast<unsigned char> (s[i]);
    if (std::isspace (c) || c == '(' || c == ')' || c == ';' || c == '"' || c == '\\') return true;
    if (i > 0 && !std::isdigit (c)) indexLike = false;
  }
  if (indexLike) return true;
  if (std::isdigit (static_cast<unsigned char> (s[0])) || s[0] == '-' || s[0] == '+') return true;
  static const char *keywords[] = { "true", "false", "and", "or", "not", "fireable", "reach", "invariant",
      "deadlock", "bound" };
  for (const char *k : keywords) if (s == k) return true;
  return false;
}

inline void printSexprName (std::ostream &os, const std::string &s)
{
  if (!needsQuotes (s)) {
    os << s;
    return;
  }
  os << '"';
  for (char c : s) {
    if (c == '"' || c == '\\') os << '\\';
    os << c;
  }
  os << '"';
}

inline void printSexprPlace (std::ostream &os, size_t p, const std::vector<std::string> *pnames)
{
  if (pnames) printSexprName (os, (*pnames)[p]);
  else os << "p" << p;
}

/** The weighted sum of places of an atom: a place, (* k place), or (+ ...). */
inline void printSexprForm (std::ostream &os, const LinearAtom &a, const std::vector<std::string> *pnames)
{
  if (a.terms.empty ()) {
    os << "0";
    return;
  }
  if (a.terms.size () > 1) os << "(+ ";
  for (size_t i = 0; i < a.terms.size (); ++i) {
    if (i > 0) os << " ";
    long long c = a.terms[i].second;
    if (c != 1) os << "(* " << c << " ";
    printSexprPlace (os, a.terms[i].first, pnames);
    if (c != 1) os << ")";
  }
  if (a.terms.size () > 1) os << ")";
}

inline void printSexpr (std::ostream &os, const LinearAtom &a, const std::vector<std::string> *pnames)
{
  os << "(" << to_string (a.op) << " ";
  printSexprForm (os, a, pnames);
  os << " " << a.constant << ")";
}

inline void printSexpr (std::ostream &os, const Expression &e, const std::vector<std::string> *pnames)
{
  switch (e.kind) {
  case Expression::Kind::True: os << "true"; break;
  case Expression::Kind::False: os << "false"; break;
  case Expression::Kind::Not:
    os << "(not ";
    printSexpr (os, e.children[0], pnames);
    os << ")";
    break;
  case Expression::Kind::And:
  case Expression::Kind::Or:
    os << (e.kind == Expression::Kind::And ? "(and" : "(or");
    for (const auto &c : e.children) {
      os << " ";
      printSexpr (os, c, pnames);
    }
    os << ")";
    break;
  case Expression::Kind::Atom: printSexpr (os, e.atom, pnames); break;
  }
}

/** One form per property; an unsupported property becomes a comment line. */
inline void printSexpr (std::ostream &os, const Property &p, const std::vector<std::string> *pnames)
{
  switch (p.kind) {
  case PropertyKind::Reachability:
  case PropertyKind::Invariant:
    os << (p.kind == PropertyKind::Reachability ? "(reach " : "(invariant ");
    printSexprName (os, p.name);
    os << " ";
    printSexpr (os, p.body, pnames);
    os << ")";
    break;
  case PropertyKind::Deadlock:
    os << "(deadlock ";
    printSexprName (os, p.name);
    os << ")";
    break;
  case PropertyKind::Bound:
    os << "(bound ";
    printSexprName (os, p.name);
    os << " ";
    printSexprForm (os, p.boundForm (), pnames);
    if (p.boundHint >= 0) os << " " << p.boundHint;
    os << ")";
    break;
  case PropertyKind::Unsupported:
    os << "; " << p.name << " unsupported: " << p.comment;
    break;
  }
}

} // namespace petri::expr

#endif /* PETRI_EXPR_SEXPRPRINTER_H_ */
