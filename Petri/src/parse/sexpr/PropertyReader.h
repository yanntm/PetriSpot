/*
 * PropertyReader.h
 *
 * S-expression forms to expr::Property values over a net. Grammar in
 * README.md. Comparisons are linearised into one LinearAtom; place and
 * transition references are indices (p<i>, t<i>) or names. A (ctl NAME f)
 * form reads a CTL formula: booleans and comparisons as elsewhere, deadlock,
 * and the path operators (EX f) (AX f) (EF f) (AF f) (EG f) (AG f)
 * (EU f g) (AU f g) (EW f g) (AW f g), case-insensitive.
 */
#ifndef PETRI_PARSE_SEXPR_PROPERTYREADER_H_
#define PETRI_PARSE_SEXPR_PROPERTYREADER_H_

#include <cctype>
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "core/Log.h"
#include "core/SparsePetriNet.h"
#include "expr/CtlFormula.h"
#include "expr/Property.h"
#include "parse/NetResolver.h"
#include "parse/sexpr/Sexpr.h"

namespace petri::sexpr
{

template<typename T>
  class PropertyReader
  {
    using Expression = petri::expr::Expression;
    using LinearAtom = petri::expr::LinearAtom;
    using Property = petri::expr::Property;
    using PropertyKind = petri::expr::PropertyKind;
    using Cmp = petri::expr::Cmp;
    using CtlFormula = petri::expr::CtlFormula;
    using CtlOp = petri::expr::CtlOp;

    /** An integer expression: sum of coeff * place, plus a constant. */
    struct LinearForm
    {
      std::vector<std::pair<size_t, long long>> terms;
      long long constant = 0;
    };

    const NetResolver<T> &net;

    [[noreturn]] static void fail (const Datum &d, const std::string &what)
    {
      throw std::string ("line " + std::to_string (d.line ()) + ": " + what);
    }

    static bool isInteger (const std::string &s)
    {
      size_t i = (s.size () > 1 && (s[0] == '-' || s[0] == '+')) ? 1 : 0;
      if (i >= s.size ()) return false;
      for (; i < s.size (); ++i)
        if (!std::isdigit (static_cast<unsigned char> (s[i]))) return false;
      return true;
    }

    /** "p12" / "t12" style index; nullopt when the atom is not of that shape. */
    static std::optional<size_t> indexRef (const Datum &d, char prefix)
    {
      const std::string &s = d.text ();
      if (d.isQuoted () || s.size () < 2 || s[0] != prefix) return std::nullopt;
      for (size_t i = 1; i < s.size (); ++i)
        if (!std::isdigit (static_cast<unsigned char> (s[i]))) return std::nullopt;
      return std::stoull (s.substr (1));
    }

    size_t place (const Datum &d) const
    {
      if (!d.isAtom ()) fail (d, "place reference expected");
      if (auto idx = indexRef (d, 'p')) {
        if (*idx >= net.placeCount ()) fail (d, "place index out of range: " + d.text ());
        return *idx;
      }
      if (auto p = net.place (d.unquoted ())) return *p;
      fail (d, "unknown place: " + d.text ());
    }

    size_t transition (const Datum &d) const
    {
      if (!d.isAtom ()) fail (d, "transition reference expected");
      if (auto idx = indexRef (d, 't')) {
        if (*idx >= net.transitionCount ()) fail (d, "transition index out of range: " + d.text ());
        return *idx;
      }
      if (auto t = net.transition (d.unquoted ())) return *t;
      fail (d, "unknown transition: " + d.text ());
    }

    static void addScaled (LinearForm &acc, const LinearForm &f, long long k)
    {
      for (const auto &t : f.terms) acc.terms.emplace_back (t.first, t.second * k);
      acc.constant += f.constant * k;
    }

    LinearForm readInt (const Datum &d) const
    {
      LinearForm f;
      if (d.isAtom ()) {
        if (!d.isQuoted () && isInteger (d.text ())) {
          f.constant = std::stoll (d.text ());
        } else {
          f.terms.emplace_back (place (d), 1);
        }
        return f;
      }
      const std::string &h = d.head ();
      const auto &items = d.items ();
      if (h == "+") {
        for (size_t i = 1; i < items.size (); ++i) addScaled (f, readInt (items[i]), 1);
      } else if (h == "-") {
        if (items.size () < 2) fail (d, "(- ...) needs an operand");
        if (items.size () == 2) {
          addScaled (f, readInt (items[1]), -1);
        } else {
          addScaled (f, readInt (items[1]), 1);
          for (size_t i = 2; i < items.size (); ++i) addScaled (f, readInt (items[i]), -1);
        }
      } else if (h == "*") {
        if (items.size () != 3) fail (d, "(* k e) takes two operands");
        LinearForm a = readInt (items[1]);
        LinearForm b = readInt (items[2]);
        if (!a.terms.empty () && !b.terms.empty ()) fail (d, "product of two place expressions is not linear");
        if (a.terms.empty ()) addScaled (f, b, a.constant);
        else addScaled (f, a, b.constant);
      } else {
        fail (d, "integer expression expected, got " + (h.empty () ? std::string ("a list") : h));
      }
      return f;
    }

    static bool cmpOf (const std::string &h, Cmp &out)
    {
      if (h == "==" || h == "=") out = Cmp::EQ;
      else if (h == "!=") out = Cmp::NEQ;
      else if (h == "<=") out = Cmp::LE;
      else if (h == ">=") out = Cmp::GE;
      else if (h == "<") out = Cmp::LT;
      else if (h == ">") out = Cmp::GT;
      else return false;
      return true;
    }

    Expression readBool (const Datum &d) const
    {
      if (d.isAtom ()) {
        if (d.text () == "true") return Expression::constant (true);
        if (d.text () == "false") return Expression::constant (false);
        fail (d, "boolean expected, got " + d.text ());
      }
      const std::string &h = d.head ();
      const auto &items = d.items ();
      Cmp op;
      if (h == "and" || h == "or") {
        std::vector<Expression> kids;
        for (size_t i = 1; i < items.size (); ++i) kids.push_back (readBool (items[i]));
        return h == "and" ? Expression::makeAnd (std::move (kids)) : Expression::makeOr (std::move (kids));
      } else if (h == "not") {
        if (items.size () != 2) fail (d, "(not e) takes one operand");
        return Expression::makeNot (readBool (items[1]));
      } else if (h == "fireable") {
        std::vector<size_t> ts;
        for (size_t i = 1; i < items.size (); ++i) ts.push_back (transition (items[i]));
        return net.anyFireable (ts);
      } else if (cmpOf (h, op)) {
        if (items.size () != 3) fail (d, "comparison takes two operands");
        LinearForm l = readInt (items[1]);
        LinearForm r = readInt (items[2]);
        LinearAtom a;
        for (const auto &t : l.terms) a.addTerm (t.first, t.second);
        for (const auto &t : r.terms) a.addTerm (t.first, -t.second);
        a.constant = r.constant - l.constant;
        a.op = op;
        a.normalize ();
        return Expression::makeAtom (std::move (a));
      }
      fail (d, "boolean expression expected, got " + (h.empty () ? std::string ("a list") : h));
    }

    static bool temporalOf (std::string h, CtlOp &op, size_t &arity)
    {
      for (auto &c : h) c = static_cast<char> (std::toupper (static_cast<unsigned char> (c)));
      arity = 1;
      if (h == "EX") op = CtlOp::EX;
      else if (h == "AX") op = CtlOp::AX;
      else if (h == "EF") op = CtlOp::EF;
      else if (h == "AF") op = CtlOp::AF;
      else if (h == "EG") op = CtlOp::EG;
      else if (h == "AG") op = CtlOp::AG;
      else {
        arity = 2;
        if (h == "EU") op = CtlOp::EU;
        else if (h == "AU") op = CtlOp::AU;
        else if (h == "EW") op = CtlOp::EW;
        else if (h == "AW") op = CtlOp::AW;
        else return false;
      }
      return true;
    }

    static bool allPreds (const std::vector<CtlFormula> &fs)
    {
      for (const auto &f : fs) if (f.op != CtlOp::Pred) return false;
      return true;
    }

    /** A CTL formula; a subtree without temporal operator or deadlock is one predicate. */
    CtlFormula readCtl (const Datum &d) const
    {
      if (d.isAtom ()) {
        if (d.text () == "deadlock") return CtlFormula::leaf (CtlOp::Deadlock);
        return CtlFormula::predicate (readBool (d));
      }
      const std::string &h = d.head ();
      const auto &items = d.items ();
      CtlOp op;
      size_t arity;
      if (temporalOf (h, op, arity)) {
        if (items.size () != arity + 1) fail (d, "(" + h + " ...) takes " + std::to_string (arity) + " operand(s)");
        if (arity == 1) return CtlFormula::unary (op, readCtl (items[1]));
        return CtlFormula::binary (op, readCtl (items[1]), readCtl (items[2]));
      }
      if (h == "and" || h == "or" || h == "not") {
        std::vector<CtlFormula> kids;
        for (size_t i = 1; i < items.size (); ++i) kids.push_back (readCtl (items[i]));
        if (h == "not" && kids.size () != 1) fail (d, "(not e) takes one operand");
        if (allPreds (kids)) {
          std::vector<Expression> es;
          for (auto &k : kids) es.push_back (std::move (k.pred));
          if (h == "not") return CtlFormula::predicate (Expression::makeNot (std::move (es[0])));
          return CtlFormula::predicate (h == "and" ? Expression::makeAnd (std::move (es)) : Expression::makeOr (std::move (es)));
        }
        if (h == "not") return CtlFormula::unary (CtlOp::Not, std::move (kids[0]));
        return CtlFormula::nary (h == "and" ? CtlOp::And : CtlOp::Or, std::move (kids));
      }
      return CtlFormula::predicate (readBool (d));
    }

    Property readForm (const Datum &d) const
    {
      const std::string &h = d.head ();
      const auto &items = d.items ();
      if (!d.isList () || items.size () < 2 || !items[1].isAtom ()) {
        fail (d, "expected (reach NAME e), (invariant NAME e), (deadlock NAME), (bound NAME e [k]) or (ctl NAME f)");
      }
      Property p;
      p.name = items[1].unquoted ();
      if (h == "deadlock") {
        if (items.size () != 2) fail (d, "(deadlock NAME) takes no body");
        p.kind = PropertyKind::Deadlock;
      } else if (h == "bound") {
        if (items.size () != 3 && items.size () != 4) fail (d, "(bound NAME e [k]) takes a form and an optional bound");
        LinearForm f = readInt (items[2]);
        if (f.constant != 0) fail (d, "a bound form is a weighted sum of places, without constant");
        LinearAtom a;
        for (const auto &t : f.terms) a.addTerm (t.first, t.second);
        a.normalize ();
        a.op = Cmp::GE;
        p.kind = PropertyKind::Bound;
        p.body = Expression::makeAtom (std::move (a));
        if (items.size () == 4) {
          if (!items[3].isAtom () || !isInteger (items[3].text ()) || std::stoll (items[3].text ()) < 0)
            fail (d, "the known bound must be a non-negative integer");
          p.boundHint = std::stoll (items[3].text ());
        }
      } else if (h == "reach" || h == "invariant") {
        if (items.size () != 3) fail (d, "(" + h + " NAME e) takes one body");
        p.kind = h == "reach" ? PropertyKind::Reachability : PropertyKind::Invariant;
        p.body = readBool (items[2]);
      } else if (h == "ctl") {
        if (items.size () != 3) fail (d, "(ctl NAME f) takes one formula");
        p.kind = PropertyKind::CTL;
        p.ctl = readCtl (items[2]);
      } else {
        fail (d, "unknown property form: " + h);
      }
      return p;
    }

  public:
    explicit PropertyReader (const NetResolver<T> &resolver)
        : net (resolver)
    {
    }

    std::vector<Property> read (const std::vector<Datum> &forms) const
    {
      std::vector<Property> props;
      props.reserve (forms.size ());
      for (const Datum &d : forms) props.push_back (readForm (d));
      return props;
    }

    std::vector<Property> read (std::string_view text) const
    {
      try {
        return read (parse (text));
      } catch (const ParseError &e) {
        throw std::string (e.what ());
      }
    }
  };

/** Read a property file. Throws std::string on a syntax or resolution error. */
template<typename T>
  std::vector<petri::expr::Property> loadProperties (const std::string &filename,
                                                     const SparsePetriNet<T> &net)
  {
    auto time = std::chrono::steady_clock::now ();
    petri::writeToLog ("Parsing property file : " + filename);
    std::ifstream in (filename, std::ios::binary);
    if (!in) throw std::string ("Cannot open property file: " + filename);
    std::stringstream buf;
    buf << in.rdbuf ();
    std::string text = buf.str ();
    NetResolver<T> resolver (net);
    PropertyReader<T> reader (resolver);
    std::vector<petri::expr::Property> result = reader.read (text);
    petri::writeToLog (
        "Parsed " + std::to_string (result.size ()) + " properties in "
            + std::to_string (std::chrono::duration_cast<std::chrono::milliseconds> (
                std::chrono::steady_clock::now () - time).count ()) + " ms.");
    return result;
  }

} // namespace petri::sexpr

#endif /* PETRI_PARSE_SEXPR_PROPERTYREADER_H_ */
