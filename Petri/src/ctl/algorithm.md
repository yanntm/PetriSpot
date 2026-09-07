# `ctl/` — the search

## Configurations

A configuration is `(s, f)`: a marking and a node of the formula in negation
normal form (`expr/CtlSimplify.h`: `Not` only on leaves, the operators
`EX AX EF AF EG AG EU AU EW AW`). Its verdict is `1`, `0` or unknown, and the
memo keeps every decided configuration and, for the unknown ones, the budget
at which they were tried, so a round with a larger budget tries them again.

## solve(s, f)

* leaf (`Pred`, `Deadlock`, `NoDeadlock`): evaluated on the marking and the
  enabled set. A formula without temporal operator is a leaf whatever its
  shape (the simplifier merges booleans over leaves into one predicate).
* `And` / `Or`: state children first, then the temporal ones; `0` / `1`
  short-circuit, an unknown child leaves the node unknown unless another
  child decides it.
* `EX g`: the enabled transitions in random order, each fired on a copy of
  the cursor, `solve(s', g)` until a `1`; all `0` is `0`; a deadlock has no
  successor and is `0`.
* `AX g`: every successor must be `1`; one `0` is `0` with that successor
  as counter-example; a deadlock is `1`.
* `E[a U b]`, `E[a W b]` (`EF b` is `E[true U b]`, `EG a` is `E[a W false]`):
  a hunt (`Hunt.h`). If it closes, `1` with the path. Else the *dual*
  `A[not b W (not a and not b)]` (`A[not b U ...]` for `W`) is tried as a
  region (`Region.h`): if the region proves it, `0`. Else unknown.
* `A[a U b]`, `A[a W b]` (`AF b`, `AG a`): the cheap refutation first, a hunt
  for the dual `E[not b U (not a and not b)]` (`E[not b W ...]` for `U`): a
  witness is `0`. Then the region; if it closes, `1`; if it meets a
  counter-example (a state where `a` and `b` both fail, a deadlock or a cycle
  where `b` never holds, for `U`), `0` with the path. Else unknown. The hunt
  comes first because an A node is mostly met inside a hunt for its parent,
  at a state where it usually fails, and a region that fails to close costs
  its whole budget.

The dual of a node is computed once and cached; a dual is never dualised
again (the flag `allowDual`), so the two methods run once each per node.

## Propagation

As in TAPAAL's dependency graphs, a verdict found for one configuration
decides others for free, and the memo receives them (`Checker::remember`):
a witness path of an E node proves the node at every state along it (the
suffix is a witness); a closed region proves the A node at every state of the
region; a counter-example path of an A node refutes it at every state along
it. A hunt from a state the path of an earlier hunt went through is thus
answered by the memo at its first step, which is what keeps `AG (EF p)`
affordable: the region enumerates the states, and `EF p` is a hunt only from
the states no earlier path crossed.

## Hunt (E nodes)

From `s`, runs of at most `runLength` steps, restarted from `s`, until the
step budget is spent. At each state: `solve(s', b)`, and `1` closes the hunt
with the path so far; else `solve(s', a)` must be `1` or the run ends. For
`W`, a deadlock closes the hunt (`a` holds there), and so does a marking seen
earlier in the run (a lasso: the loop index is recorded). The choice of the
transition is uniform, or with probability `1 - epsilon` the best of a few
sampled successors by the marking distance (`expr/Distance.h`) to `suf(b)`,
the state predicate implying `b` when the formula has one. A `0` at the start
state (`a` and `b` both fail, or `b` fails at a deadlock for `U`) is returned
as `0` at once.

## Region (A nodes)

A DFS from `s` on one cursor, firing and reverting in place, with the enabled
list copied per frame. A state where `solve(b)` is `1` is a leaf. Else
`solve(a)` must be `1`: a `0` is a counter-example (the DFS stack is the
path), an unknown aborts the region as unknown. For `U`, a deadlock or a
back edge to a state on the stack is a counter-example. A visited state that
is off the stack is not re-expanded. The region aborts as unknown when it
holds more than the state budget.

## Budgets by depth

A hunt or a region at the root of the search gets the round's budget; one
opened from a state met along a hunt or in a region (a nested obligation)
gets a hundredth of it, floored at a hundred steps and a hundred states. A state is thus probed cheaply, most probes fail fast, and the
rounds grow the probe tenfold with the rest: on CloudOpsManagement the stuck
region behind a non-live transition is proved with a few states once the
hunt lands there, while a full budget per probe spent the whole clock on
the states where the transition was still reachable.

## Rounds and reporting

The driver (`cli/CtlDriver.h`) simplifies each property, evaluates it at the
initial marking, and runs `solve` in rounds with the hunt steps and the
region states multiplied by ten per round (a thousand each in the first, six
rounds by default), under a wall clock per property. The rounds start small
on purpose: a conjunction of many obligations (the Liveness examination as
one formula) is probed child by child at the first round's price before
any child gets a large budget, the `and` stopping at the first `0` and
skipping the unknown children until the next round.
A decided property prints `FORMULA name TRUE|FALSE TECHNIQUES EXPLICIT
CTL_WALK`; with `--trace` the evidence tree follows, with transition names.
Unknown properties print an `UNKNOWN` line with `--printUnknown`.

Not in this version: threads, the verifier of the evidence tree (paths are
exact since the cursor is exact, and lassos compare exact markings; the tree
is printed as built), the LP from a state, the portfolio strategies.
