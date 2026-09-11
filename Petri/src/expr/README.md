# `expr/` — property AST

* `algorithm.md` — shape of the tree, the normal form, property kinds.
* `Expression.h` — `LinearAtom`, `Expression` (plain data), evaluation over a
  marking, printing in infix form.
* `Simplify.h` — `simplify(e)`: negation pushing, constant folding,
  flattening; the normal form described in `algorithm.md`.
* `Property.h` — `Property`, `PropertyKind`, `goal`, verdict helpers; a
  `CTL` property carries its `CtlFormula`.
* `CtlFormula.h` — the CTL tree: state predicates and the deadlock atom at
  the leaves, booleans, `EX AX EF AF EG AG EU AU EW AW`; evaluation of state
  formulas, printing.
* `CtlSimplify.h` — negation normal form (`ctlNnf`), the simplifier with the
  collapsing rules of nested operators (`ctlSimplify`, `ctlNormalize`), the
  negation of a normal-form formula (`ctlDual`), and the state predicates
  `ctlNow(f)` (implied by `f`) and `ctlSuf(f)` (implying `f`) that steer a
  walk. The checker is `ctl/`.
* `InitialState.h` — what the initial marking alone decides: the
  three-valued truth of a CTL formula there (`initialValue`), the rewrite of
  an until whose left side fails initially (`initialRewrite`), the
  requalification of `EF p` / `AG p` / `EF deadlock` into the reachability
  kinds (`requalify`), and `simplifyInitial(net, properties)` running all of
  it over a property list; called before any engine, and again after a
  reduction changed the net.
* `Hint.h` — `ParikhHint`, side information attached to a property.
* `SexprPrinter.h` — properties and expressions in the s-expression syntax of
  `parse/sexpr/` (indices, or names quoted when needed).
* `Distance.h` — TAPAAL-style estimated distance of a marking to an
  expression (zero iff it holds).

Markings are natural numbers and the normal form relies on it: leading
coefficients are positive and atoms decided by non-negativity are folded (see
`algorithm.md`). Anything consuming atoms downstream (distances, relaxed plan
goal places) may assume that form.
