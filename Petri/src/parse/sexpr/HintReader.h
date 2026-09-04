/*
 * HintReader.h
 *
 * The --hints file: s-expression forms (parikh NAME (t k)...) naming
 * properties, read into ParikhHint values and attached to the properties of
 * the same name.
 */
#ifndef PETRI_PARSE_SEXPR_HINTREADER_H_
#define PETRI_PARSE_SEXPR_HINTREADER_H_

#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "core/Log.h"
#include "core/SparsePetriNet.h"
#include "expr/Hint.h"
#include "expr/Property.h"
#include "parse/NetResolver.h"
#include "parse/sexpr/Sexpr.h"

namespace petri::sexpr
{

template<typename T>
  class HintReader
  {
    const NetResolver<T> &net;

    [[noreturn]] static void fail (const Datum &d, const std::string &what)
    {
      throw std::string ("line " + std::to_string (d.line ()) + ": " + what);
    }

    static bool isInteger (const std::string &s)
    {
      if (s.empty ()) return false;
      for (char c : s)
        if (!std::isdigit (static_cast<unsigned char> (c))) return false;
      return true;
    }

    size_t transition (const Datum &d) const
    {
      if (!d.isAtom ()) fail (d, "transition reference expected");
      const std::string &s = d.text ();
      if (!d.isQuoted () && s.size () > 1 && s[0] == 't' && isInteger (s.substr (1))) {
        size_t idx = std::stoull (s.substr (1));
        if (idx >= net.transitionCount ()) fail (d, "transition index out of range: " + s);
        return idx;
      }
      if (auto t = net.transition (d.unquoted ())) return *t;
      fail (d, "unknown transition: " + s);
    }

    /** The (TREF INT) pairs of a form, from its third item on. */
    std::vector<std::pair<size_t, long long>> pairs (const Datum &d) const
    {
      std::vector<std::pair<size_t, long long>> out;
      const auto &items = d.items ();
      for (size_t i = 2; i < items.size (); ++i) {
        const Datum &p = items[i];
        if (!p.isList () || p.items ().size () != 2 || !p.items ()[1].isAtom () || !isInteger (p.items ()[1].text ()))
          fail (p, "expected (transition non-negative-integer)");
        out.emplace_back (transition (p.items ()[0]), std::stoll (p.items ()[1].text ()));
      }
      return out;
    }

  public:
    explicit HintReader (const NetResolver<T> &resolver)
        : net (resolver)
    {
    }

    std::unordered_map<std::string, petri::expr::ParikhHint> read (const std::vector<Datum> &forms) const
    {
      std::unordered_map<std::string, petri::expr::ParikhHint> hints;
      for (const Datum &d : forms) {
        const std::string &h = d.head ();
        if (!d.isList () || d.items ().size () < 3 || !d.items ()[1].isAtom ())
          fail (d, "expected (parikh NAME (t k)...)");
        if (h != "parikh") fail (d, "unknown hint form: " + h);
        hints[d.items ()[1].unquoted ()].counts = pairs (d);
      }
      return hints;
    }
  };

/** Read a hints file; throws std::string on error. */
template<typename T>
  std::unordered_map<std::string, petri::expr::ParikhHint> loadHints (const std::string &filename,
                                                                     const SparsePetriNet<T> &net)
  {
    petri::writeToLog ("Parsing hints file : " + filename);
    std::ifstream in (filename, std::ios::binary);
    if (!in) throw std::string ("Cannot open hints file: " + filename);
    std::stringstream buf;
    buf << in.rdbuf ();
    NetResolver<T> resolver (net);
    try {
      return HintReader<T> (resolver).read (parse (buf.str ()));
    } catch (const ParseError &e) {
      throw std::string (e.what ());
    }
  }

/** Attach hints to the properties of the same name; returns how many were attached. */
inline size_t attachHints (std::vector<petri::expr::Property> &props,
                           const std::unordered_map<std::string, petri::expr::ParikhHint> &hints)
{
  size_t attached = 0;
  std::unordered_map<std::string, bool> used;
  for (auto &p : props) {
    auto it = hints.find (p.name);
    if (it == hints.end ()) continue;
    p.hint = it->second;
    used[p.name] = true;
    ++attached;
  }
  for (const auto &kv : hints)
    if (!used.count (kv.first)) std::cerr << "Warning: hint for unknown property " << kv.first << std::endl;
  return attached;
}

} // namespace petri::sexpr

#endif /* PETRI_PARSE_SEXPR_HINTREADER_H_ */
