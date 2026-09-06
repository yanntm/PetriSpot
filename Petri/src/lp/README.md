# `lp/` — in-house linear programming over the state equation

Header-only, templated on the integer type `T` of the net like the rest of
the tree. `algorithm.md` describes the method; this file maps the sources.

* `LpProblem.h` — the statement of a linear program: columns with bounds,
  ranged rows of integer coefficients (`lo ≤ a·x ≤ hi`, `±∞` allowed), an
  optional objective. A builder used by `StateEquation`; no arithmetic.
* `Simplex.h` — the solver: bounded-variable primal simplex in doubles,
  two phases, dense product-form basis inverse, partial pricing, Harris
  ratio test, pivot and time limits. Returns `LpResult` (status, solution,
  the infeasible row). Knows nothing of Petri nets.
* `StateEquation.h` — from a net and a goal expression to `LpProblem`s: the
  marking as affine forms of the transition counts, non-negativity rows,
  the atoms of one conjunct as rows, the goal's disjunctive normal form as a
  list of problems, tightening of strict comparisons.
* `Parikh.h` — a solution into `expr::ParikhHint`: rounding and the repair
  of the marking's non-negativity.
* `Refiner.h` — the `Accept / Cut / Split` interface between a candidate
  solution and the refinements (traps, read arcs, predecessor, integrality:
  planned, see `algorithm.md` section 4) and the loop that applies them.

The driver is `cli/LpDriver.h` (`--lp`, `--lpHints=FILE`): one solve per
property, `FORMULA ... TECHNIQUES STATE_EQUATION` on infeasibility, the
Parikh vectors as `(parikh NAME (t k)...)` forms in the hints file that
`--hints` reads back.

Conventions: every file under 500 lines, one responsibility each; a solver
instance and a problem belong to one thread, the net is shared read-only.
