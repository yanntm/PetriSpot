# `expr/` — property AST

* `algorithm.md` — shape of the tree, the normal form, property kinds.
* `Expression.h` — `LinearAtom`, `Expression` (plain data), evaluation over a
  marking, printing in infix form.
* `Simplify.h` — `simplify(e)`: negation pushing, constant folding,
  flattening; the normal form described in `algorithm.md`.
* `Property.h` — `Property`, `PropertyKind`, `goal`, verdict helpers.
