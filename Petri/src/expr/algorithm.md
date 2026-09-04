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
* `Unsupported`: parsed but out of the fragment (nested temporal operators,
  bounds, LTL, CTL); `comment` says why.

`goal(property)` gives the state predicate the engine has to reach, and
`verdictIfReached(property)` the MCC answer (TRUE for reachability, FALSE for
invariants).
