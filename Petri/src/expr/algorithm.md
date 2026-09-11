# `expr/` — property AST: shape, normal form, simplification

## Shape

A property body is a boolean tree over **linear atoms**. The AST is plain
data (no virtual dispatch, no template metaprogramming); evaluation is a small
template function over any marking type exposing `get(place) -> integer`.

```
Expression ::= True | False
             | Not(Expression)
             | And(Expression+) | Or(Expression+)
             | Atom(LinearAtom)
LinearAtom ::= sum_i coeff_i * m[place_i]   cmp   constant      cmp ∈ {==, !=, <=, >=, <, >}
```

Terms of an atom are sorted by place, coefficients are non-zero integers,
each place appears once. The constant and the accumulated value are 64-bit.

Nothing in the AST refers to a concrete syntax. Fireability (`is-fireable(t)`
in MCC) is desugared by the parser into `And(m[p] >= w for (p,w) in pre(t))`
so it never reaches the AST. Deadlock is not an atom either: it is a property
kind handled by the exploration engine.

## Normal form (`simplify`)

`simplify(e)` returns an equivalent expression such that:

1. no `Not` node remains: negation is pushed to atoms (`cmp` is negated) and
   through `And`/`Or` by De Morgan;
2. `<` and `>` are rewritten to `<=` and `>=` on integers (`x < k` is
   `x <= k-1`);
3. the leading coefficient of an atom is positive (the whole atom is negated
   otherwise: `-p <= -1` becomes `p >= 1`); an atom with no term is folded to
   `True`/`False`, and so is an atom whose coefficients are all positive when
   the non-negativity of markings decides it (`p + q <= -1` is `False`,
   `p >= 0` is `True`);
4. `And`/`Or` are flattened (a child of the same kind is spliced in), neutral
   children are dropped, an absorbing child collapses the node, duplicate
   children are removed, a single child replaces the node, an empty `And` is
   `True` and an empty `Or` is `False`.

The simplifier is total and idempotent. It does not attempt implication
reasoning between atoms (`x <= 3 && x <= 5`), which belongs to a later,
optional pass.

## Properties

A `Property` has a name, a kind and a body:

* `Reachability`: body is `phi` from `EF phi`; a state satisfying `phi` makes
  the property TRUE.
* `Invariant`: body is `psi` from `AG psi`; a state violating `psi` makes the
  property FALSE. The exploration goal is `simplify(Not(psi))`.
* `Deadlock`: `EF deadlock`; body unused.
* `Bound`: maximise a weighted sum of places (the terms of the single atom in
  `body`); `boundHint` is a known upper bound or -1. The exploration goal is
  `form >= hint` when there is a hint, and there is no reachable goal without
  one: the engine reports the largest value seen.
* `CTL`: a formula with nested path operators, kept as a `CtlFormula` in
  `ctl`; checked by `ctl/` after normalisation (`CtlSimplify.h`).
* `Unsupported`: parsed but out of the fragment (LTL, unknown elements);
  `comment` says why.

## Initial state (`InitialState.h`)

Before any exploration, every property is confronted with the initial marking
`M0`, after Bonneland et al. (Petri Nets 2018, section 3) as ITS-Tools applies
it. A reachability body true at `M0` is TRUE, an invariant body false at `M0`
is FALSE, the deadlock property is TRUE when nothing is enabled at `M0`. A CTL
formula, in normal form, is given a three-valued truth at `M0`:

* a state formula is evaluated (`deadlock` from the enabled transitions);
* `EG f`, `AG f` are false when `f` is false at `M0`; `EF f`, `AF f` are true
  when `f` is true there; `EX`, `AX` say nothing;
* an until or weak until `[a U b]` is true when `b` holds at `M0`, false when
  neither `a` nor `b` holds there.

A decided formula becomes its constant. Otherwise the formula is rewritten
with what `M0` decides: the children of a boolean are rewritten in turn, and
an until whose left side fails at `M0` is replaced by its right side, which
then has to hold at `M0` itself. Operands under a path operator speak of other
states and are left alone. The result is simplified again and *requalified*:
`EF p` with `p` a state predicate becomes a `Reachability` property, `AG p` an
`Invariant`, `EF deadlock` the `Deadlock` property. So `E[a U AG p]` with `a`
false initially is answered as the invariant `AG p`, by the reachability
engines and under the reachability reductions, not by the CTL checker.

The pass runs on the parsed properties before the structural reduction reads
their kinds and supports, and again on the reduced net after its constant
places were substituted (`reduction/PropertyFacts.h`).

## CTL normal form (`ctlNormalize`)

Negation is pushed to the leaves (`not E[a U b] = A[not b W (not a and not
b)]`, `not EX = AX not`, `not EF = AG not`, and so on), predicates are
simplified as above and merged when they are siblings under a boolean, and
nested operators collapse: `EF EF a`, `EF AF a`, `AF EF a`, `EF E[a U b]`,
`EF A[a U b]`, `A[a U EF b]`, `E[a U EF b]` are all `EF` of the innermost
operand; `AF AF a = AF a`, `A[a U AF b] = AF b`, `EG EG a = EG a`,
`EG AG a = AG a`, `AG AG a = AG a`; `AG (a and b) = AG a and AG b`,
`EF (a or b) = EF a or EF b`, an `EF` disjunct is pulled out of the right
side of an until or an `AF`; `deadlock` and `not deadlock` as an until's
operand reduce it; `X` of a constant becomes the deadlock atom (`EX true`
is `not deadlock`, `AX false` is `deadlock`: a deadlock has no successor).
The rules are TAPAAL's (verifypn, `CTL-formula-equivalence-rewriting.pdf`).

`goal(property)` gives the state predicate the engine has to reach, and
`verdictIfReached(property)` the MCC answer (TRUE for reachability, FALSE for
invariants).
