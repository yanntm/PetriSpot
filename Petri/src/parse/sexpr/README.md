# `parse/sexpr/` — properties as s-expressions

The tool-to-tool property syntax (see `INTEROP.md` section 4): one form per
property, indices or names for places and transitions, nothing to escape but
quotes. The reader is the libHSC `datum` reader (`hsc::surface::sexpr`),
copied so the two projects keep one syntax without depending on each other.

* `Sexpr.h` — `Datum` (an atom or a list, with its source line), `parse(text)`
  to a vector of top-level forms, `write(os, datum)`. Atoms are runs of
  non-delimiter characters or `"quoted strings"` (kept with their quotes in
  the atom text, `\"` and `\\` escapes); `;` comments to end of line.
* `PropertyReader.h` — forms to `expr::Property` values over a net:
  `(reach NAME BEXP)`, `(invariant NAME BEXP)`, `(deadlock NAME)`. Booleans
  are `true`, `false`, `(and ...)`, `(or ...)`, `(not e)`, `(CMP e e)` with
  `CMP` in `== != <= >= < >`, `(fireable t...)`; integers are literals, place
  references, `(+ ...)`, `(- a b...)` or `(- a)`, `(* k e)`. A comparison is
  linearised into one `LinearAtom`; a product of two non-constant sides is an
  error. `p<digits>` and `t<digits>` are indices, `"quoted"` atoms are names,
  any other atom is a name. `loadProperties(file, net)` reads a file.

Errors are thrown as `std::string` with the line of the offending form, like
the MCC loader. The printer for this syntax is `expr/SexprPrinter.h`, and the
choice between MCC XML and s-expressions is `parse/PropertyFile.h`.
