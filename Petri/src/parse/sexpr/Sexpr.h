/*
 * Sexpr.h
 *
 * S-expression reader: text to a tree of Datum (atom or list, with source
 * line). Adapted from libHSC's hsc::surface::sexpr so both projects read the
 * same syntax; quoted-string atoms added. Knows nothing of properties.
 */
#ifndef PETRI_PARSE_SEXPR_SEXPR_H_
#define PETRI_PARSE_SEXPR_SEXPR_H_

#include <cctype>
#include <cstddef>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace petri::sexpr
{

/** One node: an atom (symbol, integer or quoted string token) or a list. */
class Datum
{
  bool list_ = false;
  std::string text_;
  std::vector<Datum> items_;
  int line_ = 0;

public:
  static Datum atom (std::string text, int line)
  {
    Datum d;
    d.text_ = std::move (text);
    d.line_ = line;
    return d;
  }
  static Datum list (std::vector<Datum> items, int line)
  {
    Datum d;
    d.list_ = true;
    d.items_ = std::move (items);
    d.line_ = line;
    return d;
  }

  bool isAtom () const
  {
    return !list_;
  }
  bool isList () const
  {
    return list_;
  }
  int line () const
  {
    return line_;
  }
  /** Token text of an atom (quotes included for a quoted string); empty for a list. */
  const std::string& text () const
  {
    return text_;
  }
  /** Children of a list; empty for an atom. */
  const std::vector<Datum>& items () const
  {
    return items_;
  }
  /** Head symbol of a list "(head ...)", or "" when empty or not a list. */
  const std::string& head () const
  {
    static const std::string empty;
    if (!list_ || items_.empty () || !items_.front ().isAtom ()) return empty;
    return items_.front ().text_;
  }
  bool isQuoted () const
  {
    return !list_ && text_.size () >= 2 && text_.front () == '"';
  }
  /** Inner text of a quoted atom with escapes resolved; the text itself otherwise. */
  std::string unquoted () const
  {
    if (!isQuoted ()) return text_;
    std::string out;
    for (size_t i = 1; i + 1 < text_.size (); ++i) {
      if (text_[i] == '\\' && i + 2 < text_.size ()) ++i;
      out += text_[i];
    }
    return out;
  }
};

class ParseError : public std::runtime_error
{
  int line_;
public:
  ParseError (int line, const std::string &what)
      : std::runtime_error ("line " + std::to_string (line) + ": " + what), line_ (line)
  {
  }
  int line () const
  {
    return line_;
  }
};

namespace detail
{

/** Cursor over the source, tracking the line for diagnostics. */
class Scanner
{
  std::string_view src;
  size_t pos = 0;
  int line_ = 1;

  static bool isDelimiter (char c)
  {
    return std::isspace (static_cast<unsigned char> (c)) || c == '(' || c == ')' || c == ';' || c == '"';
  }

public:
  explicit Scanner (std::string_view s)
      : src (s)
  {
  }

  /** Skip whitespace and ';' comments; false at end of input. */
  bool skipTrivia ()
  {
    while (pos < src.size ()) {
      char c = src[pos];
      if (c == '\n') {
        ++line_;
        ++pos;
      } else if (std::isspace (static_cast<unsigned char> (c))) {
        ++pos;
      } else if (c == ';') {
        while (pos < src.size () && src[pos] != '\n') ++pos;
      } else {
        return true;
      }
    }
    return false;
  }
  char peek () const
  {
    return src[pos];
  }
  int line () const
  {
    return line_;
  }
  void advance ()
  {
    ++pos;
  }

  /** A bare atom: up to the next delimiter. */
  std::string atomToken ()
  {
    size_t start = pos;
    while (pos < src.size () && !isDelimiter (src[pos])) ++pos;
    return std::string (src.substr (start, pos - start));
  }

  /** A quoted atom, quotes and escapes kept verbatim. */
  std::string quotedToken ()
  {
    int startLine = line_;
    size_t start = pos;
    ++pos; // opening quote
    while (pos < src.size ()) {
      char c = src[pos];
      if (c == '\\' && pos + 1 < src.size ()) {
        pos += 2;
        continue;
      }
      if (c == '\n') ++line_;
      ++pos;
      if (c == '"') return std::string (src.substr (start, pos - start));
    }
    throw ParseError (startLine, "unterminated string");
  }
};

inline Datum readDatum (Scanner &s)
{
  int line = s.line ();
  if (s.peek () == '(') {
    s.advance ();
    std::vector<Datum> items;
    for (;;) {
      if (!s.skipTrivia ()) throw ParseError (line, "unterminated list");
      if (s.peek () == ')') {
        s.advance ();
        return Datum::list (std::move (items), line);
      }
      items.push_back (readDatum (s));
    }
  }
  if (s.peek () == ')') throw ParseError (line, "unexpected ')'");
  if (s.peek () == '"') return Datum::atom (s.quotedToken (), line);
  return Datum::atom (s.atomToken (), line);
}

} // namespace detail

/** Parse the top-level forms of text. Throws ParseError. */
inline std::vector<Datum> parse (std::string_view text)
{
  detail::Scanner s (text);
  std::vector<Datum> forms;
  while (s.skipTrivia ()) forms.push_back (detail::readDatum (s));
  return forms;
}

/** Render one datum, no newline; write(parse(t)) round-trips modulo whitespace. */
inline void write (std::ostream &os, const Datum &d)
{
  if (d.isAtom ()) {
    os << d.text ();
    return;
  }
  os << '(';
  bool first = true;
  for (const Datum &item : d.items ()) {
    if (!first) os << ' ';
    write (os, item);
    first = false;
  }
  os << ')';
}

} // namespace petri::sexpr

#endif /* PETRI_PARSE_SEXPR_SEXPR_H_ */
