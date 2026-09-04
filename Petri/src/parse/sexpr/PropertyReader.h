/*
 * PropertyReader.h
 *
 * S-expression forms to expr::Property values over a net. Grammar in
 * README.md. Comparisons are linearised into one LinearAtom; place and
 * transition references are indices (p<i>, t<i>) or names.
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

    Property readForm (const Datum &d) const
    {
      const std::string &h = d.head ();
      const auto &items = d.items ();
      if (!d.isList () || items.size () < 2 || !items[1].isAtom ()) {
        fail (d, "expected (reach NAME e), (invariant NAME e) or (deadlock NAME)");
      }
      Property p;
      p.name = items[1].unquoted ();
      if (h == "deadlock") {
        if (items.size () != 2) fail (d, "(deadlock NAME) takes no body");
        p.kind = PropertyKind::Deadlock;
      } else if (h == "reach" || h == "invariant") {
        if (items.size () != 3) fail (d, "(" + h + " NAME e) takes one body");
        p.kind = h == "reach" ? PropertyKind::Reachability : PropertyKind::Invariant;
        p.body = readBool (items[2]);
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
