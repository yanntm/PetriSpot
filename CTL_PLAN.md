# An explicit CTL checker on the walkers: witnesses first (design, 2026-09-07, to comment)

The question: can PetriSpot's explicit machinery (sparse walkers, portfolio,
quests, LP over the state equation, invariants) settle CTL formulas whose
answer has a *small witness*, one way or the other, and leave exhaustive
proofs to the symbolic engine (`its-ctl` on libDDD, the only CTL backend of
ITS-Tools; libHSC has no CTL either)? This note says what a CTL witness is,
what TAPAAL's explicit engine does that we should borrow, what exists here,
what is missing, and proposes a folder `ctl/` with a lean first milestone:
basic walks that try to close a formula by a witness. Nothing is built yet.

## 1. What a CTL witness is

MCC formulas are evaluated at the initial marking over the reachability graph.
With negation pushed to the atoms (CTL has a negation normal form once the
weak until `W` is admitted) a formula is a tree of

```
phi ::= atom | phi and phi | phi or phi
      | EX phi | AX phi
      | E[phi U phi] | A[phi U phi]      (EF b = E[true U b], AF b = A[true U b])
      | E[phi W phi] | A[phi W phi]      (EG a = E[a W false], AG a = A[a W false])
```

with `not E[a U b] = A[not b W (not a and not b)]` and `not A[a U b] =
E[not b W (not a and not b)]`. Exactly one of `phi`, `not phi` is led by an
existential operator; the engine searches both and reports whichever closes
(a witness of `not phi` is the counter-example of `phi`).

Evidence for a formula at a state `s` is a finite tree:

* `atom`: the marking; `and` / `or`: evidence for all / one child;
* `EX phi`: one enabled transition and evidence for `phi` at its successor;
* `E[a U b]`: a firing sequence from `s`, evidence for `a` at every state
  before the last, for `b` at the last;
* `E[a W b]`: the same, or a sequence ending on a *lasso* (a state already on
  the sequence) or on a *deadlock*, with `a` all along;
* `AX phi`: evidence for `phi` at every successor: the enabled list, cheap;
* `A[a U b]`, `A[a W b]`: the *region* of states reachable from `s` through
  states where `b` fails, with evidence for `a` at each; for `U` the region is
  finite, acyclic and has no deadlock. Exhaustive, so cheap only when the
  region is small or a structural argument closes it without enumeration.

"Strong at counter-examples, not at proofs" therefore means: E nodes are
searched by walkers, A nodes are closed by cheap proofs (structural, or an
enumeration under a budget), and the engine answers only when the whole tree
closes. Anything else is UNKNOWN, for the symbolic engine.

Semantics assumed, to confirm against the MCC manual before coding: a
deadlock loops on itself, so `EG a` holds at a deadlock satisfying `a`,
`AX phi` at a deadlock is `phi` there, and `A[a U b]` fails at a deadlock
where `b` fails. Deadlocks are witness-friendly and the walkers find them.

## 2. The shape of MCC formulas (measured)

The 16 formulas of `CTLCardinality` and `CTLFireability` on the models under
`bench/models/` (a throwaway script; a kept one is phase 0): 3 to 5 path
quantifiers deep, 4 to 15 quantifiers per formula, every operator present, E
and A interleaved (AirplaneLD-PT-0010 CTLFireability: `EF AG EG AF AU`,
`AG AN EG AF EU EF`, one `EF EG`, one `AF AG`). Plain `EF` / `AG` are a
handful. So:

* a checker for the top-level fragment answers almost nothing new; the tree
  is searched as a whole;
* A nodes sit *under* E paths (`E[a U AG c]` needs an `AG` proof at the end
  of a hunted path) and E nodes sit under A regions (`AG (EF c)` needs a hunt
  from every region state): both nestings must be cheap where they can be and
  budgeted where they cannot;
* on the small models of the campaign `its-ctl` already solves 285 of 288
  CTLC and 273 of 288 CTLF (HANDOFF 2026-09-07). The explicit engine's value
  is on the models where the symbolic one blows up, the regime where the
  walkers beat the reachability engines. Expect a *fraction* of the formulas
  there, and measure that fraction first (phase 0).

## 3. What TAPAAL does, and what we take from it

TAPAAL's `verifypn` checks CTL explicitly as a local fixed point over a
*dependency graph* (Liu and Smolka's linear local algorithm, 1998; the
extended graphs with negation edges and the *certain-zero* optimisation of
Dalsgaard, Enevoldsen, Fogh, Jensen, Jensen, Johannsen, Larsen, Muñiz, Olesen
and Srba, Petri Nets 2016 and 2017; the same engine is distributed over
threads in their later work). The ideas, stated in our words:

* A **configuration** is a subgoal `(s, node)`: a marking and a subformula.
  Its value is 1 (holds), 0 (fails) or unknown.
* A **hyperedge** from a configuration to a set of configurations is one way
  to prove it: `(s, a and b) -> {(s,a),(s,b)}`; `(s, a or b)` has two
  hyperedges of one target each; `(s, EX phi)` has one hyperedge per enabled
  transition, `(s, AX phi)` one hyperedge with every successor;
  `(s, E[a U b])` has `{(s,b)}` and, per successor `s'`, `{(s,a), (s',E[a U b])}`;
  `(s, A[a U b])` has `{(s,b)}` and one hyperedge `{(s,a)} ∪ {(s',A[a U b])
  for all s'}`. `EG`/`AG` go through negation edges.
* The **local algorithm** expands configurations on demand from the root,
  keeps a waiting list, and propagates 1 backwards the moment a hyperedge
  has all its targets at 1. The configurations touched by that propagation
  *are the witness tree*: the algorithm is a witness search, and it never
  enumerates more of the graph than the witness needs (plus what it explored
  in vain).
* **Certain zero** propagates definite failure: a configuration is certainly
  0 when every hyperedge has a target that is certainly 0. For an A node one
  failing successor suffices: the counter-example side, again without
  exhaustive exploration. For an E node all hyperedges must fail, which is
  the exhaustive side.
* The **search order** of the waiting list is a heuristic (they use the
  distance of the marking to the atoms, the same estimate as their
  reachability engine): which configuration to expand next is what decides
  whether the witness is found in thousands of steps or never.

What we take: configurations as subgoals, the memo of values, 1 and
certain-zero propagation, and the observation that the witness is the closed
fragment of the graph. What we change:

* An E node is **never expanded into all its hyperedges**. A walker picks one
  successor at a time, pursuing `(s', E[a U b])` as a *quest*: a chain of
  hyperedges of the same node along the walk, with `(s, a)` checked at each
  step. The graph never learns that an E node fails, which is exactly the
  exhaustive question we do not ask; it only learns that it holds. TAPAAL
  expands hyperedges on the waiting list; we expand along a walk, thousands
  of configurations a millisecond, keeping only what the witness needs.
* An A node's hyperedge **is expanded in full, under a budget**, or closed by
  a structural proof (LP over the state equation from `s`, section 4) that
  TAPAAL does not have.
* The search order is the portfolio's: strategies, quests, shares, restarts,
  the pool, over many open configurations at once as `TargetSet` does today
  over many targets.
* No partial order reduction: PetriSpot has none, and LoLA's CTL-preserving
  stubborn sets are their ground. If regions, not hunts, turn out to be the
  bottleneck, that is the time to look.

The quest framework already has this shape: `QuestSweep` picks a target,
runs `sync` toward it from where the processes stand, and picks again from
the state where the quest ended; `ComponentStrategy` stages subgoals
(barriers, their pre-places) on a stack and pops them as they close. A CTL
formula is that stack given by the formula rather than by the net: reach a
state promising `b`, from there close the sub-obligation, and if it does not
close, walk on.

## 4. What exists here and plugs in

* **Hunting a state predicate** from a marking: `Portfolio` over a
  `TargetSet`, every strategy, hints, the pool, the scheduler. An E node whose
  right side is a state predicate is today's reachability target.
* **Deadlocks**: `--findDeadlock`, `DeadlockDistance`, saturation. An `EG a`
  ending on a deadlock is a deadlock hunt under a constraint.
* **Successors**: `EnabledSet` gives the enabled list, `Marking::peek`
  applies and reverts an effect. `AX` and `EX` cost `|enabled|` peeks.
* **Global truths about atoms**: the invariant engine (semiflows in-process,
  as `sync` does) and the LP over the state equation (`lp/`: `Simplex.h`,
  `StateEquation.h`). An atom whose violation is infeasible under
  `m = m0 + C x, x >= 0` holds on every reachable marking: `AG atom` is TRUE
  everywhere, `EF not atom` FALSE everywhere, and the CTL simplifier folds
  the node. The same LP *from a state `s`* (`m = s + C x`) proves `AG atom`
  at `s` with no enumeration: the cheapest universal leaf we can have, and
  the one the on-the-fly checkers do not lean on.
* **Traces and verification**: recorded traces replayed by an independent
  walker before printing. A witness tree is a tree of such traces.
* **Time sharing**: `Scheduler`, `Coordinator`, resumable tasks. Every piece
  of the search below is a resumable task.

## 5. What is missing

1. **The formula.** `expr/` stops at state predicates; the MCC parser tags a
   nested quantifier Unsupported. Needed: a CTL AST in NNF over
   `expr::Expression` leaves, the MCC XML and s-expression parsers extended,
   a CTL simplifier (booleans, `X` on constants, atoms decided globally by the
   LP folded, `E[a U b]` with `b` true at `s`), and the two approximations of
   a node by a state predicate that steer the walkers: `suf(node)`, a
   predicate implying it (`suf(atom) = atom`, `suf(E[x U y]) = suf(y)`,
   `suf(A[x U y]) = suf(y)`, `suf(E[x W y]) = suf(y)`, `suf(and)` / `suf(or)`
   pointwise, `false` for `X`, `AG`, `EG`), and `now(node)`, a predicate it
   implies (`now(E[x U y]) = now(x) or now(y)`, `now(EG x) = now(x)`,
   `true` for `X`).
2. **Path constraints.** `E[a U b]` needs a walk that stays in `a`: a guard in
   the walker refusing a transition whose successor violates `a`. Sparse:
   only transitions touching a place of `a` in the direction that can break
   it need the peek, the `up`/`down` lists of `TargetIndex` encode this.
3. **Cycle detection** for lassos of `E[a W b]`: a rolling hash of the
   marking (`sum z[p] * m[p]` mod 2^64, updated per touched place,
   O(|effect|)) and a run-local set of hashes with their step index, cleared
   at reset. A match is a lasso *candidate*; the verifier confirms it on exact
   markings when replaying, so a collision costs a replay, never a verdict.
   Eight bytes per step of the run.
4. **Refusable claims.** When `b` is temporal, reaching `now(b)` opens a
   sub-obligation `(s', b)` rather than closing the target; if it fails the
   target stays open and `s'` is remembered as tried. `TargetSet::claim` is a
   CAS today; it becomes a callback that may refuse.
5. **Budgeted region exploration** for `A[a U b]`, `A[a W b]`: a DFS from `s`
   through `not b` states with apply/revert (allocation free, sparse), a
   visited set of hashed markings, a cap on states; back edges and deadlocks
   are counter-examples for `U`, allowed for `W`; `(s', a)` is a subgoal at
   each state. Resumable: its state is the DFS stack.
6. **The configuration memo** `(marking, node) -> 1 | 0 | unknown(budget)`,
   shared by the subgoals: TAPAAL's assignment, restated with budgets.
   Without it `AG (EF c)` re-hunts `c` from every region state.
7. **Witness trees** and their verifier: a tree of (firing sequence, loop
   index or deadlock, children), checked by a plain recursive evaluator over
   enumerated markings before any `FORMULA` line.

## 6. The search

`solve(s, node, budget) -> 1 | 0 | unknown`, memoised:

* **atom**: evaluate (globally decided atoms are constants already);
* **and / or**: children in order of estimated cost (atoms, `X`, hunts,
  regions); an `or` runs its hunts as one target set;
* **EX phi**: enabled transitions in heuristic order (the successor's
  distance to `now(phi)`), `solve(s', phi)` until 1;
* **AX phi**: every enabled transition, all must be 1; a deadlock is
  `solve(s, phi)`;
* **E[a U b]**, **E[a W b]**: a hunt. The walker runs from `s` under guard
  `a`, targets `suf(b)` when it is not `false`, else `now(b)`; on arrival it
  opens `solve(s', b)` and a refusal keeps the target open; for `W` it also
  claims on a deadlock and a lasso candidate. A hunt is a portfolio task with
  a step budget; a Parikh hint to `suf(b)` from `s` is today's `--hints`
  mechanism with `s` in place of `m0`;
* **A[a U b]**, **A[a W b]**: proofs in order of cost: `b` at `s`; for `W`
  with `a` a state predicate, the LP from `s` proving `AG a`, or a deadlock
  at `s` satisfying `a`; then the region DFS under the budget. Over budget is
  unknown for this `(s, node)` at this budget and the enclosing hunt walks
  on: it now looks for *another* state where the proof closes. That is the
  heuristic content of the A side: steer the walker to states where the
  universal obligation is cheap, and the deadlock distance (few enabled
  transitions, a terminal region nearby) is a first estimate of cheapness.

Budgets escalate by rounds (steps and region states tenfold per round, as
`WalkDriver` does under `--escalate`), the memo keeping what closed at lower
budgets; `phi` and `not phi` are solved in alternation on the same memo. The
witness tree is assembled from the memo's evidence when the root closes,
verified, and printed: the MCC `FORMULA` line, and with `--trace` the tree as
a nested list of firing sequences with loop marks, which is also what the
tests check.

Threads: hunts are portfolio tasks already; region DFSs and `AX` fan-outs
become tasks of the same scheduler; the memo is the shared structure (a
striped hash map, or per thread with merge; decide when measuring). The first
milestone is single-threaded.

## 7. Integration

A driver mode on the same binary (`--ctl`, or the property kind decides), fed
the MCC XML or the s-expression syntax. On the Java side it runs *beside*
`its-ctl` as `ParallelWalk` runs the reachability walk beside the diagram
engines: whoever closes a formula publishes into `DoneProperties`, and the
explicit engine's UNKNOWN costs the portfolio nothing. The Java presolving
(`AtomicReducer`, knowledge facts, reductions) stays upstream and hands the
engine a reduced net and a simplified formula, as for LTL.

## 8. Risks, honestly

* **Yield.** The fraction of MCC CTL formulas with a small witness tree is
  unknown; phase 0 measures it. If most A nodes under E paths need regions of
  millions of states the engine adds little to what reductions and the
  reachability fragment already give. The models where the symbolic engine
  fails are the ones that count; count on those.
* **Soundness.** A wrong MCC verdict is expensive and the portfolio takes a
  verdict as final. Only trees that pass the independent verifier are
  printed; region proofs are re-run by the verifier; hashes propose, never
  decide.
* **Memory.** Visited sets and the memo are bounded by the budgets; markings
  are stored sparse; the memo evicts unknown entries by size.
* **Semantics.** Self-loop on deadlocks, `X` at a deadlock, `is-fireable` in
  a guard: confirm against the MCC manual and `its-ctl` verdicts before the
  first release; a differential test against `its-ctl` on the bench models is
  the acceptance test.
* **Depth.** Hunts under regions under hunts multiply budgets; the rounds and
  the memo are the defence, a per-formula wall clock the last resort.

## 9. Plan of attack

0. **Measure before building** (half a day). A script in `Petri/test/`
   classifies the formulas of the bench models by shape after NNF: E-led or
   A-led, depth, A nodes under E paths and the nature of their right sides,
   and, once the simplifier exists, what remains after folding globally
   decided atoms. Set against the `its-ctl` verdicts and times of the
   campaign it says where the explicit engine could matter, and whether to
   go on.
1. **`ctl/` formula** (about 600 lines): AST in NNF, MCC XML and s-expression
   parsing, CTL simplifier with LP-decided atoms, `now` / `suf`, printing.
   Folder `README.md` and `algorithm.md` first.
2. **Basic walks, single thread** (about 1 000 lines): `solve` with the memo,
   `EX` / `AX` fan-out, region DFS, hunts on the existing `Walker` with the
   guard, the rolling hash and refusable claims (walk-side changes, about 200
   lines), random and best-first strategies only, witness trees and the
   verifier. Validated on `Petri/examples/` and small bench nets against
   `its-ctl`. This is the milestone that tells whether the idea pays.
3. **Hunts through the portfolio**: quests toward `suf` / `now`, deadlock
   distance for the `W` / `G` ends, LP hints from a state, rounds and
   escalation. Measure the yield on the bench models.
4. **LP proofs from a state** for `W` / `G` leaves, shared memo, hunts and
   regions as scheduler tasks on several threads.
5. **MCC driver and Java side**: `--ctl` end to end beside `its-ctl`,
   differential test on the bench, then the cluster on CTLCardinality and
   CTLFireability.

Out of scope: fairness (none in MCC CTL), CTL*, exhaustive CTL (that is
`its-ctl`), partial order reduction for the regions.
