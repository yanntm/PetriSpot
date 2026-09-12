# Small inequality examples

These ordinary P/T nets exercise phase-1 harvesting in the existing executable:

```
build/petri64 -i Petri/examples/inequalities/net1.pnml --Pflows --collectInequalities
```

Use `--Psemiflows` for the positive basis path, or `--Tflows` to collect
transition vectors with one-sign effects. Collection is incomplete and follows
the chosen elimination path. These models do not request extra candidate search.

| Net | Transitions | Initial marking |
| --- | --- | --- |
| 1 | p0 -> p1; p1 -> p0 + p2 | p0 = 1 |
| 2 | p0 -> p1; p1 + p2 -> p0 | p0 = 1, p2 = 5 (N) |
| 3 | p0 -> p1; p1 -> nothing | p0 = 1 |
| 4 | nothing -> p0; p0 -> p1 | all zero |
| 5 | nothing -> p0; p0 -> p1; p1 -> nothing | all zero |

All arcs have weight one. Transition names are t0, t1, t2 in the listed order.
Place inequalities print their initial-marking constants. The matrix API and
KERS outputs return coefficients and direction only.

## Observed phase-1 collection

With default heuristics, the P-flow runs collect:

| Net | Collected inequalities | Equality basis |
| --- | --- | --- |
| 1 | p2 >= 0 | p0 + p1 = 1 |
| 2 | p2 <= 5; p0 + p2 <= 6 | p0 + p1 = 1 |
| 3 | p0 <= 1 | empty |
| 4 | p1 >= 0; p0 + p1 >= 0 | empty |
| 5 | none | empty |

The lower bounds starting at zero also certify nondecreasing sums, even though
their displayed constants alone follow from token nonnegativity.

P-semiflow collection differs: net 1 gives p1 + p2 >= 0; net 2 gives only
p0 + p2 <= 6; net 4 gives only p1 >= 0. Nets 3 and 5 give the same inequalities
as the P-flow runs. Existing equality bases agree on these examples.

This shows the intended incompleteness: net 3's p0 + p1 <= 1 is missed. The
singleton equation for the sink can discard p1 before a total-count candidate
is formed. We do not change pivots or retain discarded candidates to recover it.
For net 5 there is no nonzero nonnegative-weight place sum with a uniform effect
sign: the source, transfer and sink force both weights to zero for either sign.

On transitions, T-flow harvesting finds no one-sign vectors for nets 1 and 2,
despite t0 + t1 having respectively positive and negative net effect. It finds
t1 and t0 + t1 with negative effects for net 3, t0 with positive effect for
net 4, and t2 and t1 + t2 with negative effects for net 5. Net 5 additionally
has the ordinary T-semiflow t0 + t1 + t2. Again, discarded opportunities depend
on the existing elimination order.

T-semiflow harvesting gives the same transition certificates except on net 3,
where it finds only t1. Neutral vectors are never included in either inequality
matrix; they remain exclusively in the equality results.

To collect both orientations in one run:

```
build/petri64 -i Petri/examples/inequalities/net2.pnml --Pflows --Tflows --collectInequalities
```

For binary exports, use separate P and T invocations so each output matrix has
one unambiguous variable space.

## Projection consumer

The accompanying `.sexpr` files ask small bound-related queries. libHSC's
approximation consumes harvested inequalities by default, for example:

```
hsc-pn -i Petri/examples/inequalities/net2.pnml --props Petri/examples/inequalities/net2.sexpr \
  --approx 2 --approx-only --printUnknown
```

Net 2's p2 becomes covered: the resource inequalities are proved and the exact
maxima p2=5 and p0+p2=6 are reported (both attained initially). Without
harvesting, these queries remain UNKNOWN because p2 is projected away. Net 3
now proves p0<=1 and maximum p0=1, but its total bound remains uncovered. See
the consumer observations in `../../../INEQUALITIES.md`. Add
`--no-approx-inequalities` to compare against equality-only approximation.
