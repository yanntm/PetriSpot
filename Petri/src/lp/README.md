# `lp/` — in-house linear programming over the state equation

Header-only, templated on the integer type `T` of the net like the rest of
the tree. `algorithm.md` describes the method; this file maps the sources.

* `LpProblem.h` — the statement of a linear program: columns with bounds,
  ranged rows of integer coefficients (`lo ≤ a·x ≤ hi`, `±∞` allowed), an
  optional objective. A builder used by `StateEquation`; no arithmetic.
* `Basis.h` — the basis inverse in product form: a signed diagonal and one
  sparse eta per pivot; `ftran`, `btran`, `pivot`, `reset`.
* `Simplex.h` — the solver: bounded-variable primal simplex in doubles,
  slack starting basis, two phases, candidate pricing refreshed by full
  sweeps, Harris ratio test, Bland after a run of degenerate pivots, periodic
  rebuild of the basis, pivot and time limits; `solve(base, extraRows)` for a
  branch or a cut without copying the base. Returns `LpResult`. Knows nothing
  of Petri nets. `--lpDebug` traces the rebuilds with a residual check.
* `StateEquation.h` — from a net and a goal expression to `LpProblem`s: the
  marking as affine forms of the transition counts, non-negativity rows,
  the atoms of one conjunct as rows, the goal's disjunctive normal form as a
  list of problems, tightening of strict comparisons.
* `Parikh.h` — a solution into `expr::ParikhHint`: rounding and the repair
  of the marking's non-negativity.
* `Refiner.h` — the `Accept / Cut / Split` interface between a candidate
  solution and the refinements (traps, read arcs, predecessor, integrality:
  planned, see `algorithm.md` section 4) and the loop that applies them.
* `DeadlockRefiner.h` — the dead-marking condition built lazily: a split on
  the pre-places of a transition still enabled in the candidate.

The driver is `cli/LpDriver.h` (`--lp`, `--lpHints=FILE`, `--lpTime=S`,
`--lpSolves=N`): one solve per property (a refinement tree for a deadlock), `FORMULA ... TECHNIQUES STATE_EQUATION` on infeasibility, the
Parikh vectors as `(parikh NAME (t k)...)` forms in the hints file that
`--hints` reads back.

Conventions: every file under 500 lines, one responsibility each; a solver
instance and a problem belong to one thread, the net is shared read-only.
