# `ctl/` — explicit CTL checking, witnesses first

A local CTL checker over the explicit reachability graph, built on the walk
substrate (`walk/WalkNet.h`, `Marking.h`, `EnabledSet.h`). It closes a formula
when it finds a witness or a counter-example of bounded size and answers
UNKNOWN otherwise; it never proves anything by exhausting the state space.
`algorithm.md` describes the search; the design and its context are in
`CTL_PLAN.md` at the repository root. The formula AST and its simplifier live
in `expr/` (`CtlFormula.h`, `CtlSimplify.h`); the parsers in `parse/`.

* `Checker.h` — `Checker<T>`: `solve(cursor, formula)` over a memo of
  `(marking, node)` verdicts; `EX`/`AX` by successor enumeration, `E` nodes
  by hunts, `A` nodes by budgeted regions, each also trying the dual method
  on its negation; verdicts propagated along witness paths and closed
  regions (`remember`); budgets and the clock in `Verdict.h`.
* `Verdict.h` — the three-valued verdict, the budgets of a round, the
  counters printed per property.
* `Hunt.h` — the hunt for `E[a U b]` / `E[a W b]` from a state: guarded
  random walks (epsilon-greedy on the distance to `suf(b)`) with restarts,
  deadlock and lasso detection for `W`.
* `Region.h` — the region DFS for `A[a U b]` / `A[a W b]`: the states
  reachable through `not b`, all owing `a`, with a cap on their number;
  a deadlock or a cycle inside is a counter-example for `U`.
* `Cursor.h` — a marking with its enabled set, copied for sub-obligations,
  fired and reverted in place along a walk or a DFS.
* `Evidence.h` — the witness tree (paths, loops, region sizes, children) and
  its printing with transition names.

Semantics: MCC's, as `its-ctl` and TAPAAL compute it: a deadlock ends its
path, so `EG a` holds at a deadlock satisfying `a` and `A[a U b]` fails at a
deadlock where `b` fails, while `EX phi` is false there and `AX phi` true
(no successor).
