# `lp/` — the in-house linear programming engine

An over-approximation of reachability by linear arithmetic, computed in
process on the sparse substrate of `core/`, without an external solver. It
answers three questions the walker and the total examinations ask:

* **Is this goal reachable at all?** If the relaxation is infeasible the goal
  is unreachable in the net: a sound `FALSE` for a reachability query, a
  sound `TRUE` for an invariant, a dead transition for QuasiLiveness.
* **How would one get there?** A solution of the relaxation is a vector of
  firing counts, the *Parikh vector*, handed to the walker as a hint
  (`expr/Hint.h`, the `parikh` strategy).
* **How high can this sum of places go?** The maximum of a place sum under
  the relaxation is a valid upper bound; with the walker's lower bound it
  closes an UpperBounds atom.

The relaxation is the *state equation*: any marking `M` reachable from `M0`
satisfies `M = M0 + C·x` for some `x ≥ 0` in `N^T`, `C` the incidence matrix
(effects, `flowTP − flowPT`). Over the rationals it is a linear program.
Every refinement below adds linear rows or splits the problem in two, so the
engine stays a linear program throughout and the same solver serves every
stage.

## 1. The problem, and the shape of the formulation

Variables are the transitions only: `x_t ≥ 0`, one per transition. Places do
not get variables; the marking is the affine expression
`s_p(x) = M0_p + Σ_t C[p,t]·x_t`, substituted wherever an atom mentions `p`.
This keeps the column count at `|T|` and the row count at what the question
needs:

* **non-negativity of the marking**: one row `Σ_t C[p,t]·x_t ≥ −M0_p` per
  place `p` (skipped for a place no transition touches);
* **the goal's atoms**: an atom `Σ_p a_p·s_p cmp k` becomes the row
  `Σ_t (Σ_p a_p·C[p,t])·x_t cmp k − Σ_p a_p·M0_p`, a sparse combination of
  the rows of `C`;
* **an objective**: minimise `Σ_t x_t` for a short Parikh vector, maximise
  the goal's linear form for a bound, none for plain feasibility.

Rows are ranged, `lo ≤ a·x ≤ hi`, with `±∞` where a side is free; an equality
has `lo = hi`. Coefficients are kept as the integers they are (the net's
weights are integers, the atoms' coefficients are integers): the numerical
solver reads them as doubles, an exact checker can read them as they are.

Markings are integers, so a strict inequality is tightened before it enters:
`s < k` is `s ≤ k − 1`, `s > k` is `s ≥ k + 1`, `s ≠ k` is a disjunction of
the two. The relaxation over the rationals of the tightened system still
contains every integer solution; nothing is lost in soundness and the
relaxation is sharper.

## 2. Boolean structure

A goal is a tree of `And`, `Or`, `Not` over atoms (`fireable` is expanded
into `≥` atoms at parse time, `expr/`). The linear program takes a
conjunction of atoms. The goal is put in *disjunctive normal form* with
negations pushed to the atoms (`Not(a cmp k)` is `a negate(cmp) k`), each
conjunct is one linear program, and the goal is infeasible when every
conjunct is. The number of conjuncts is bounded (`MAX_CASES`); beyond it the
engine declines the goal rather than enumerate. MCC atoms are small, the
usual goal is a handful of conjuncts, the QuasiLiveness atoms are one.

## 3. The solver

A bounded-variable primal simplex on the sparse columns, in double
precision, with the tolerances and safeguards of the trade:

* **Standard form.** Ranged rows get a slack each, `lo ≤ a·x + slack ≤ hi`
  turned into `a·x − r = 0` with `lo ≤ r ≤ hi`, so the constraint matrix is
  `[A | −I]` and every row is an equality; variables carry bounds
  `[l, u]` with `u = +∞` for a free upper side.
* **Starting basis and phase 1.** The start is `x = 0`, the initial marking
  itself, so every non-negativity row is satisfied and its slack is basic;
  only a row violated at the start (an atom) gets a basic artificial, driven
  to zero by minimising the sum of their absolute values (the two-phase
  method). Artificials of satisfied rows are fixed at zero and costless: a
  basis of fixed artificials would take hundreds of zero-step pivots to
  leave, and that alone made the first prototype ten times slower. A row
  whose artificial cannot reach zero is the certificate of infeasibility.
* **Basis.** A dense inverse of the basis (rows are few on the nets that
  matter: hundreds of places against a hundred thousand transitions) updated
  by the product form at every pivot, refactorised every `REFACTOR` pivots
  by Gaussian elimination on the basis columns. Sparse LU is the next step
  when a net with tens of thousands of places asks for it; the interface
  (`Basis`) hides the representation.
* **Pricing.** A full sweep of the columns costs milliseconds on a hundred
  thousand transitions, a pivot a fraction of one, so a sweep must serve
  many pivots: it keeps the sixty-four columns of largest reduced cost, the
  following iterations choose the best among them (Dantzig) as long as one
  is eligible, then sweep again. Under Bland's rule the first eligible
  column in index order enters, the scan stopping there. Steepest edge is
  the refinement when pricing dominates again.
* **Ratio test.** Harris two-pass with a feasibility tolerance, so a tiny
  pivot never enters the basis; bound flips for boxed variables.
* **Degeneracy.** After `DEGENERATE` consecutive zero-step pivots Bland's
  rule takes over until a pivot makes progress; it guarantees termination.
  Bound perturbation, removed at the end, is the planned refinement when
  degenerate runs show up in the pivot counts.
* **Limits.** A pivot budget and a wall-clock deadline per problem; the
  result says which limit stopped it, so a caller in a total examination
  can move on and a caller with time can retry with more.

The result of a solve is one of `Optimal`, `Infeasible`, `Unbounded`,
`PivotLimit`, `TimeLimit`, with the primal solution when there is one and,
on infeasibility, the row whose artificial stayed positive.

### Soundness of what is published

A published infeasibility is a proof about the net, so a floating-point
verdict is not enough on its own. The engine keeps the integer data of the
problem, and an *exact checker* (`ExactCheck`, on `core/Rational.h` over
the integer type of the binary) recomputes the final basis in rational
arithmetic and confirms the certificate: for infeasibility a Farkas vector
`y` with `y·A ≥ 0`, `y·b < 0` (read off the phase 1 duals), for optimality
the reduced costs. Until the checker exists (it is the next piece), the
`--lp` driver prints the floating-point verdict as a `FORMULA` line for the
tests, and nothing in the walk or the total examinations consumes it. A
Parikh vector needs no check: it is a hint, the walk that follows it is the
proof.

### Parikh vector from a solution

The solution is rational; the hint wants integers. Each `x_t` is rounded to
the nearest integer (a fractional firing counts as one); rounding up where
the rounded vector no longer keeps the marking non-negative is planned. The
hint keeps the transitions with a positive count, `(parikh NAME (t k)...)`,
and the walker does the rest (`walk/ParikhStrategy.h`).

## 4. Refinement: cuts and splits

The state equation admits spurious solutions (a token borrowed before it is
produced, a trap emptied and refilled). The refiners of ITS-Tools express
these as constraints for an SMT solver; here each is either a linear row
added to the program or a split of the program in two, decided from a
candidate solution. The interface is `Refiner::examine(problem, solution)`,
returning `Accept`, `Cut(rows)` or `Split(branches)`; the driver loops until
a solution is accepted, the problem is infeasible, or the budget is spent.

* **Traps** (`TrapRefiner`, planned). Let `Q` be the places empty in the
  candidate marking. The greatest trap inside `Q` is the fixpoint that
  removes from `Q` every place some consumer of which produces nothing in
  the set. If it holds a place marked in `M0`, the candidate is spurious and
  the cut `Σ_{p∈S} s_p ≥ 1` is a linear row; a greedy minimisation of `S`
  gives the sharper cut. No solver is needed to find the trap.
* **Read arcs** (`ReadFeedRefiner`, planned). A transition `t` reading a
  place `p` below its initial marking needs a feeder to have fired before it.
  On a candidate with `x_t > 0` and too little feeding, the problem splits:
  `x_t = 0` in one branch, `Σ_{feeders} x ≥ delta` in the other.
* **Predecessor** (planned). For an invariant, the last transition of a
  path to a violating marking touched the goal's support, its post-set is
  covered and the goal held false before it: one branch per candidate last
  transition, each a few rows, worth it when the support is touched by tens
  of transitions and skipped when by thousands.
* **Integrality** (planned). A GCD test on each equality row is free; a
  bounded branch and bound on the fractional variables is the escalation,
  behind the same `Split` interface.
* **Deadlock** (`DeadlockRefiner`, built). A dead marking starves every
  transition, one disjunction per transition: far too many conjuncts for a
  normal form, so the condition is built lazily. The refiner reads the
  candidate marking, takes a transition still enabled there (the one with
  the fewest pre-places) and splits on which pre-place to starve, the
  emptiest first, one row `s_p ≤ w − 1` per branch. Every branch infeasible
  is a proof of deadlock freedom. Measured on 2026-09-06: DatabaseWithMutex
  PT-02 and PT-04 proved deadlock-free; Philosophers 5, 10 and 100 give the
  deadlock's Parikh vector (every left fork once) in `n + 1` solves. The
  larger instances hit the cost of a solve, not the branching: see below.

### What the prototype measured about scale

Every solve of a branch starts from scratch with a dense `m × m` inverse,
`m` the number of places: 15 ms a solve at 830 rows (DatabaseWithMutex
PT-10, 3 966 branches in 60 s without closing), 170 ms and 200 MB at 5 000
(Philosophers 1 000), refused above 6 000 rows. Two changes fix this, in
order: the **product form of the inverse** with sparse eta vectors instead
of the dense matrix (a basis here is the slack identity plus a few hundred
structural columns, so the work per operation is `O(k·m)` for `k` structural
basics, and the row limit goes away), and a **warm start of a split**: the
parent's optimal basis stays dual feasible when a row is added, so a branch
is a handful of dual simplex pivots rather than a phase 1. For deadlock
freedom on the wide instances the branching itself may stay too wide; the
classical route is structural, a marked trap in every siphon, and belongs
to the siphon work of WALK_PLAN.md section 10.6.

Cuts on a solved program are re-solved warm: the dual simplex is the natural
tool once rows are added, and is the first extension of `Simplex` after the
prototype; until then a cut re-solves from the current basis by phase 1.

## 5. Where it plugs in

* **Driver mode.** `petri -i model --props=file --lp [--lpHints=out.sexpr]`
  solves every property once: an infeasible goal prints a `FORMULA` verdict
  (`TECHNIQUES STATE_EQUATION`), a feasible one writes its Parikh vector to
  the hints file the walker reads with `--hints`. Bounds maximise their form
  and print `BOUND` lines. This is the test harness as raw PetriSpot.
* **Inside a walk run.** Before the sweep, the goals get an LP each under a
  per-goal budget: infeasible ones are answered and leave the sweep, the
  others carry a hint into the focused rounds. On the total examinations
  this is the proved side we could not answer.
* **Program to program.** The Java side asks through the hint protocol of
  INTEROP.md: a property file in, verdicts and hints out.

## 6. Shape of the code

```
lp/
  LpProblem.h       variables with bounds, ranged rows of integer coefficients, objective; builder
  Simplex.h         the bounded-variable primal simplex over doubles: Basis, pricing, ratio test, limits
  StateEquation.h   net + goal -> LpProblem: the marking as affine forms, atoms as rows, DNF cases
  Parikh.h          a solution -> ParikhHint (rounding, non-negativity repair)
  Refiner.h         the Accept / Cut / Split interface and the refinement loop
  TrapRefiner.h     (planned) trap fixpoint and its cut
  ExactCheck.h      (planned) rational confirmation of a certificate
```

`Simplex` knows nothing of nets, `StateEquation` nothing of pivots, the
refiners see a problem and a solution. Thread safety: the net and its
incidence matrix are shared read-only, every problem and solver instance is
owned by one thread.

## 7. Yardsticks

* Every target the walker has claimed is feasible (soundness in the cheap
  direction): the 1 000 ResIsolation targets, the Erlangen and Stigmergy
  witnesses.
* Time per problem on the wide nets: ResIsolation-PT-N10P4 (445 places,
  147 855 transitions) under a second per atom for feasibility.
* Hints that help: the three Erlangen targets and the CAN gathering target
  with an LP hint against without.
* Infeasibility found where the campaign proved nothing: the QuasiLiveness
  residue of the seven families of `TOTAL_QUERIES.md`.
