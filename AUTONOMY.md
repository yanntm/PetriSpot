# Standing alone: examinations, procedures, and the loop

A design, to comment, written 2026-09-11 evening, first comments folded in.
The subject is the next big step: PetriSpot answering an MCC examination
from raw inputs to the last FORMULA line without ITS-Tools, with the ideas
ITS-Tools carries but built as one framework rather than inherited as one
`Application.java`. Not a competitor: the same author, the version written
knowing what ten years of ITS-Tools taught, with the engineering done
properly this time. ITS-Tools keeps its specificities and strengths.
Vocabulary follows `PORTFOLIO.md` (fact, goal, cohort, strategy, profile);
what is new here is how examinations, procedures and the classification of a
property fit together, and where the symbolic engine lives.

## 1. Where we stand

The ingredients exist, each in its folder with its own design file:

| ingredient | where | state |
| --- | --- | --- |
| PNML and MCC XML parsing, s-expressions, PNET | `parse/`, `io/` | mature, vendored by libHSC |
| the property tree, simplification, initial-state rules | `expr/` | mature; `InitialState.h` requalifies EF/AG to reachability kinds |
| structural reductions, the counting record, the preparation pipeline | `reduction/` | reachability and deadlock inventory complete, STATESPACE with the record, `prepare` to a fixpoint |
| state equation, simplex, refiners, dead transitions | `lp/` | built; exact checker pending |
| explicit walks, strategies, multi-target, Parikh, bounds | `walk/` | mature, the engine ITS-Tools calls |
| the CTL checker | `ctl/` | first proof of concept, section 9 of `CTL_PLAN.md` |
| the symbolic engine | libHSC, `hsc-pn` | the `.hsc` language and its decision diagrams; `--reduce` now runs our reductions in memory |
| the harness, the oracles, the campaign pages | `~/git/MCC-drivers`, `~/git/MCC-analysis` | ours, working |

What still runs only in ITS-Tools, and why it matters:

* **Orchestration**: `Application.java` decides per examination what to
  run, in what order, with what budget. It is one long method with the
  choices of ten years inside. This is what this document replaces.
* **Coloured nets**: unfolding and the skeleton (`SparseHLPetriNet`).
* **LTL**: Spot for the automata, the stutter and knowledge tests, LTSmin
  as the product engine.
* **SMT-backed reductions**: implicit places, dead transitions by
  z3/yices, with the read-arc abstraction.
* **Decision diagrams**: the ITS engine behind `MultiOrderRunner`.

Of these the last two have native replacements in progress (the LP and
libHSC), the first is this document, and coloured nets and LTL are named at
the end as the two migrations that come after.

## 2. The architecture decision: one tree, the engine as a library

hsc-pn is the Petri net face of libHSC: it parses nets and properties with
our vendored headers, reduces with our vendored kernel, and hands a `.hsc`
model to the engine. Every part of it that knows what a Petri net is belongs
here. The decision:

* **libHSC keeps what is general**: the `.hsc` command language, the
  diagram engine, the surface calculus, its own example adapters for a
  couple of other languages. It points to PetriSpot for Petri net support.
* **PetriSpot takes the Petri net side of hsc-pn**: `to_surface`,
  `props_to_surface`, the NUPN unit tree and `decompose`, `pn_solver`,
  `pn_approx` and its pass, the StateSpace values, the pumping pair. They
  become the `symbolic/` folder here, a client of libHSC's public API.
* **The dependency runs from PetriSpot to libHSC**, optional at build time.
  `petri64 -hsc` (the name is provisional) activates the symbolic
  procedures; without libHSC the same binary runs everything else. No more
  vendoring in either direction: libHSC is a library PetriSpot links, found
  by CMake as a sibling checkout or an installed package.
* **One binary answers an examination.** ITS-Tools keeps calling it while
  it still owns the coloured and LTL flows, then stops.

Why this direction and not the other: a Petri net tool that links a
symbolic engine is ordinary; a symbolic engine that embeds a Petri net tool
and its solving loop is not. And it removes the cycle that today forces the
preparation to run inside hsc-pn.

## 3. The vocabulary

**Question.** One thing to decide, with its knowledge state: a boolean
(`unknown`, TRUE, FALSE), an interval for a bound, a value for a count.
Questions come from the property file, from the examination (the four
StateSpace values, one question per place for OneSafe), or from a procedure
that spawns subquestions. A question carries its *formula*, its *kind*, its
*support* (the places it observes), its *provenance*, and what was spent on
it. This is `PORTFOLIO.md`'s goal, renamed to avoid the clash with the
reduction goal.

**Kind.** What the formula's shape says it is, after simplification, and
nothing else. The kinds form a partial order, not a chain: what one kind
observes decides which reductions are sound for it, and CTL and LTL observe
different things. CTL sees the next step (`EX`, `AX`) and branching; LTL
sees infinite paths under an implicit universal quantifier and never the
next step of one branch, so the stutter-insensitive reductions serve LTL
and the stutter-insensitive CTL fragment, the branching-preserving ones
serve CTL, and neither set contains the other.

```
constant  <  initial-state decidable  <  invariant / reachability of a state predicate
          <  reachability with fireability atoms  <  bound  <  deadlock
          <  CTL without next (stutter-insensitive)  <  CTL
    and, apart:  LTL (stutter-insensitive by construction; its own engine, section 6)
```

Kinds only ever move *down*, by rewriting: `EF p` is a
reachability question, `AG p` an invariant, `E[p U q]` with `p` false at
the initial marking is `q` at the initial marking, an until whose right side
is constant false is constant false. `expr/InitialState.h` does these today,
by the Bonneland rules. The principle to hold to: **a procedure registers
for kinds, never for examination names.** A CTL examination whose formula is
`AG p` is an invariant question and gets the invariant procedures, in their
order, with their budgets; the CTL checker only ever sees what no simpler
kind absorbed. This is what makes the simple case the natural one rather
than an edge case: there is no CTL-specific path for `AG` because there is
no CTL question left.

**Fact.** What a procedure learned: a verdict, a witness, a bound, an
invariant, a structural fact (dead transition, constant place, bound of a
place), an exploration statistic, a profile signal. Facts are stated over
the base model by place name when they concern the net, or tagged with the
generation of the reduced net that produced them otherwise; a stale fact is
dropped, not translated (`PORTFOLIO.md`, "Model, and what a fact is about").

**Knowledge.** The set of facts of a cohort, indexed by question and by
model generation. It is the shared structure between threads; the
`DoneProperties` of ITS-Tools, made explicit. A fact narrows questions: a
verdict closes one, a bound tightens an interval, a structural fact may
rewrite a formula (a dead transition makes a fireability atom false, a
constant place substitutes a value) and so lower its kind.

**Pass.** A transformation of the model: the structural reductions for a
goal, the constant substitution, the dead-transition retirement, the
unfolding of a coloured net, the abstraction that drops unbounded places.
A pass produces a new model generation and says what it preserved: the
reduction goal names the kinds and supports it is sound for, and whether
counts survive (the record). Passes are the reduction `Coordinator`'s
business today; they stay there.

**Procedure.** A decision procedure: a worker with an applicability
predicate over a question's kind and the knowledge, a cost class, what
facts it can produce, and what model shape it needs (the reduced net, the
original, an over-approximation). Examples: initial-state evaluation,
structural deadlock deduction, the state equation (infeasible proves an
invariant, a Parikh vector hints a walk), the dead-transition tests, the
walks in their strategies, the CTL checker, the symbolic reachable set, the
symbolic count, the pumping pair, k-induction, a bound by invariants.

**Examination.** A macro goal: the MCC name, the questions it opens from
the inputs, how the reduction goal follows from the open kinds, the
schedule of procedures by cost, its budget, and the output format. It is a
configuration, not code with branches: what `Application.java` does by `if
("StateSpace".equals(examination))` becomes one table row per examination.

**Cohort, conjoined goal, subgoal.** As in `PORTFOLIO.md`: the open
questions of one examination on one instance, the model generations they
live on, and the clock they share. The cohort starts as one **conjoined
goal**: every open question attacked together, on one model reduced for
their union support, by procedures that answer many at once (the
multi-target walks, one symbolic set, one state equation with many
targets); the flood of QuasiLiveness queries is the case that makes this
mandatory. When the conjoined goal is stuck, the portfolio **isolates**: a
question, or a few sharing a support, becomes a subgoal with its own model
(reduced on its support alone, where the metrics say the reduction pays),
its own budget, and every procedure of its class, its kind degrading as it
learns. Isolation is gated by memory: reduction applicability is a
model-wide fact, so a rule that found nothing on the conjoined model is not
retried on every subgoal. A round is one pass over the schedule at a cost
level, over the conjoined goal and the live subgoals.

## 4. Examinations as configurations

One row each; the schedule lists procedure families in cost order, each
gated by the kinds still open and by the profile.

| examination | questions opened | reduction goal | schedule, cheap to dear | answer |
| --- | --- | --- | --- | --- |
| StateSpace | STATES, TRANSITIONS, MAX_TOKEN_IN_PLACE, MAX_TOKEN_PER_MARKING | STATESPACE with the record; the dead-transition tests inside | invariant bounds (a MAX_TOKEN upper bound), the symbolic count on the reduced net, the pumping pair when a place is unbounded (`+inf`) | the four `STATE_SPACE` lines, only what the record vouches for |
| ReachabilityCardinality, ReachabilityFireability | one boolean per formula, kinds reachability or invariant | REACHABILITY on the union support, per-question when the profile allows | constants and initial state; structural deductions; the state equation (proofs, Parikh hints); walks by strategy under clock budgets; the symbolic set; the CTL checker as a last engine | FORMULA lines with techniques |
| ReachabilityDeadlock | one boolean | DEADLOCK | initial state; structural deduction (a source transition, a stuck net); the state equation with the deadlock refiner; walks; the symbolic set with the deadlock atom | one FORMULA line |
| UpperBounds | one interval per expression | REACHABILITY on the expression's support, bound dominance | initial value; invariants and the state equation maximum (upper); walks (lower); the symbolic maximum (exact); close when lo = hi | one value per formula |
| OneSafe, StableMarking, QuasiLiveness, Liveness | generated: one per place or transition (`TOTAL_QUERIES.md`) | REACHABILITY or LIVENESS; the dead-transition tests answer QuasiLiveness's negatives directly | the total-examination multi-target walks; structural facts as answers (a dead transition is not quasi-live, a constant place is stable); invariant bounds for OneSafe; the symbolic set for the rest | the aggregated verdict and the per-object vectors |
| CTLCardinality, CTLFireability | one per formula, **kind from the shape**: most are reachability or invariant after the initial-state rules, the rest CTL | by the open kinds: REACHABILITY when none is CTL, SI_CTL when none uses next, LTL-safe otherwise | everything the lower kinds get, then the CTL checker and the symbolic CTL for the CTL kind | FORMULA lines |
| LTLCardinality, LTLFireability | one per formula, kind LTL | SI_LTL or LTL | stays in ITS-Tools until Spot is a dependency here; then the stutter and knowledge tests become procedures over facts | FORMULA lines |

Flags configure a row rather than fork the code: `--totalTime`, a budget
per procedure family, `--no-<procedure>` and `--only-<procedure>`, `-hsc`
to admit the symbolic procedures, `--reduce` off for the degraded mode. A
named profile is a saved set of such flags; the reference configuration
reproduces today's behaviour, so an experiment is a named deviation.

## 5. The loop

The refinement loop of ITS-Tools' `ReachabilitySolver`, with its two
monotone quantities made explicit: **cost rises, observation shrinks.**

```
round 0   parse; open the questions; classify; constants; initial state; answer what is decided
round k   reduce the model for the open kinds and the union support (record when counting)
          cheap procedures on the new generation: structural deductions, dead transitions,
              the state equation; answer, rewrite, reclassify
          if the cohort is small, per-question models on each question's own support
          mid procedures under clock budgets: walks by strategy, the symbolic set with a cap
          reclassify: a closed question leaves the support; a bound narrows; an until rewrites
          if anything changed since round k began: k+1, budgets escalate
          else the dear procedures on what is left: the CTL checker, the symbolic engine
              with the rest of the clock, k-induction
stop      no open question, or the clock
```

The portfolio's job is the alternation: an hour of budget, time sliced
between the conjoined goal and the subgoals, between cheap and dear
procedures, between explicit, linear and symbolic engines, each share
following what it earned in the previous round. That alternation is where
the state of the art is beaten, not in any one engine.

### Epochs

The unit of the loop is the **epoch**: a periodic rendez-vous where the
cohort is reconsidered as a whole. At an epoch boundary:

* the goal set is re-read: what closed, what a fact from one goal says
  about another (a bound closes an invariant, a dead transition kills a
  fireability atom), what reclassified;
* finished goals return their unused resources to the pool;
* procedures are re-budgeted: the successful ones escalate, the ones that
  earned nothing are discarded for this instance, and some signals admit
  the exotic ones (a shape that failed twice admits a different order, a
  stalled walk admits k-induction);
* when few goals are left and their characteristics make them disjoint,
  support or nature, the cohort attacks them individually; when they
  batch, they batch by kind, the stutter-insensitive CTL formulas together
  and the rest together, each batch on the reductions its kind allows.

"Stuck" is relative and uses every signal, not a counter. For the symbolic
engine it usually means memory growth or states accreting far too slowly:
the shape gets a negative signal, and the failed run is not waste but an
under-approximation to harvest, reachable states with a good heuristic
score handed to the walks, or facts the set already proves. The epoch is
the brain of the operation: today a policy of rules over the profile and
the knowledge; later, when the tensors of knowledge and outcomes exist
across campaigns, a learned policy choosing the next epoch's strategies
and budgets. The design only has to make the inputs of that policy
explicit and recorded, which is what the profile is. ITS-Tools does this
implicitly at best, with hard-coded procedures and timeouts and little
knowledge shared; making it explicit is the advantage we hold.

Scheduling is **cooperative**, on our own thread pool with its park and
resume (`walk/`'s tasks and coordinator): procedures poll their deadlines,
and what enters the pool is decided centrally at the epoch. A literal
rendez-vous is the first version because it keeps the reasoning simple;
the same reasoning can later run as one more task scanning a snapshot of
the pool's load and dependencies, without stopping anyone.

What is new against the Java loop:

* **Reclassification is a step**, run after every batch of facts, and it
  feeds the reduction goal of the next round. A CTL formula that becomes an
  invariant after a rewrite changes the goal from LTL-safe to REACHABILITY,
  which admits the agglomerations, which is why the order matters.
* **Budgets come from the clock and escalate per round** (`PORTFOLIO.md`
  G2). A procedure that ran out of budget is *resumable*: it keeps a cursor
  and continues next round rather than restarting (`lp/DeadTransitions.h`
  does this already).
* **Every procedure reports its metrics into the profile**: found, tested,
  time, budget hit. The profile gates the next round: a dead-transition
  density of zero after a full sweep turns the test off for this instance,
  a walk with no new verdict for two rounds cedes its share.
* **The record travels with the model.** A counting examination keeps its
  record through the passes that maintain it and refuses the others; the
  same code answers a StateSpace and a Reachability cohort, the goal decides
  which rules run.
* **Threads share knowledge, not models.** A procedure runs on its own copy
  of a generation; facts go to the knowledge structure; the loop reads it
  between steps. The pool of `PORTFOLIO.md` is the scheduler.
* **A procedure's model comes with a set of initial states**, the initial
  marking being the singleton case. An explicit engine takes a marking, a
  symbolic one a set; a subgoal spawned from inside a search (section 8)
  hands over where it stands. This is the one generalisation the
  interfaces of section 7 carry from the start.

## 6. `Application.java`, mapped

| Java | here | state |
| --- | --- | --- |
| `MccTranslator` (parse, unfold, `rebuildSpecification`) | `parse/`; the coloured unfolding to come | PT complete, COL later |
| `ReachabilitySolver.checkInInitial` | `expr/InitialState.h` as the initial-state procedure | done |
| `ReachabilitySolver.applyReductions(rt, doSMT)` | the reduce pass with the goal from the open kinds; the SMT rules become the LP and symbolic procedures | structural done; implicit places by LP to do |
| `ReachabilitySolver` loop, `Effort` | section 5 | to build |
| `DoneProperties` | knowledge | to build, small |
| `GlobalPropertySolver` | the generated questions of the total examinations | the walks exist (`TOTAL_QUERIES.md`), the generator to move |
| `UpperBoundsSolver.applyReductions` | bound questions with intervals, bounds dominance rule | rule exists unscheduled |
| `LTLPropertySolver`, stutter and knowledge tests | procedures over facts, on our own LTL engine: the canonical syntactic omega-semigroup of the formula as a computed transition matrix, a few accepting pairs and entry points per letter, in place of a Büchi automaton (the local `LTLToBuchi` work) | later, and better than Spot's path |
| `MultiOrderRunner.runMultiITS`, `startHsc` | the symbolic procedures in `symbolic/`, in process | hsc-pn `--reduce` is the prototype |
| the examination `if` chain | the table of section 4 | to build |

Migration order, each step retiring one Java flow: StateSpace (the
symbolic count is ours already), then the reachability family with
UpperBounds, then the total examinations, then CTL. The **unfolder** comes
right after the procedure catalogue exists: porting it accurately is about
one session, the skeleton-based decisions another, and the coloured
strategies are then written on top of the catalogue rather than beside it.
**LTL** comes last and on its own engine, not on Spot: the semigroup
construction above is why it waits, and why it will not be a port. **SMT**:
the LP is the native, dependency-free prototype and it stands; when the
refiners want a real solver, one that people have spent years on is
plugged behind the same `lp/` interface rather than grown here.

## 7. Code organisation

New folders, each with its `README.md` and `algorithm.md` first:

* `exam/`: the examination table, one file per examination, and the
  factory from the MCC name and the flags. No solving code.
* `proc/`: the procedure interface and the adapters that wrap what exists
  (`InitialState`, `StateEquation`, `DeadTransitions`, the walks, the CTL
  checker, the symbolic ones) into it. Thin: the algorithms stay in their
  folders.
* `loop/`: questions, kinds and the reclassification, knowledge, budgets,
  the rounds of section 5. This is where `reduction/Pipeline.h`'s
  `prepare` dissolves into steps: classify, constants, initial state,
  reduce. The kernel below it does not move.
* `symbolic/`: the Petri net side of hsc-pn, a client of libHSC.

The interfaces, kept to what the loop needs:

```cpp
enum class Kind { Constant, Initial, Invariant, Reachability, Fireability, Bound, Deadlock, CtlSI, Ctl, Ltl };
struct Question { std::string name; expr::Property formula; Kind kind; Support support; Knowledge state; Provenance from; Spent spent; };
struct Procedure {
  virtual std::string_view name() const = 0;
  virtual Cost cost() const = 0;                                   // cheap, mid, dear
  virtual bool applies(const Question&, const Knowledge&) const = 0;
  virtual Outcome run(Cohort&, const Model&, Budget) = 0;          // facts in, resumable
};
struct Examination {
  std::vector<Question> open(const Inputs&) const;
  reduction::Goal goal(const Cohort&) const;                       // from the open kinds
  std::vector<Procedure*> schedule(Cost level) const;
  void answer(const Question&, std::ostream&) const;               // the line protocol
};
```

Not a pattern catalogue: one virtual interface for procedures because they
are many and grow, plain data for examinations and questions because they
are tables. The trace and metrics hooks of the reduction kernel are the
model for a procedure's reporting.

## 8. What "a CTL AG gets the love it deserves" means concretely

Two different things, and the design keeps them apart.

At the top, classification: today the CTL examination runs `prepare`, which requalifies `AG p` to an
invariant kind, then hands the remaining properties to the CTL driver,
which decides by kind whether to walk or to check. That is the edge-case
shape: the kind is known, the driver is CTL's. In the design above the
question enters the cohort as an invariant, the examination's schedule is
the CTL row's, but every procedure gates on the kind, so the invariant
procedures fire and the CTL checker's predicate is false for it. Nothing
in the CTL checker knows about `AG`; nothing in the reachability procedures
knows they serve a CTL examination. The same holds downward: a reachability
question the walks and the state equation left open can be handed to the
symbolic set or the CTL checker, since their predicates accept the lower
kinds too, at their cost level.

Inside a search, subgoals: a CTL checker deciding `E[p U AG q]` meets `AG q`
at some state, or in the symbolic engine at a set of states, and must
handle it internally: an invariant question from that state or set. Whether
the checker spawns such subgoals into the cohort, to be attacked by the
invariant procedures under their own budget with the set as initial states
(section 5), or keeps them as internal searches, is **open in this
architecture**; the interfaces allow both, since a procedure takes a set of
initial states, and that is all the provision made now. What is already
clear: an explicit engine given a set picks a state and walks it, a partial
answer that is often enough; whether re-reducing for a new initial state
pays depends on the model, useless on a live net, worth it only where
behaviour stabilises. A symbolic run spawning goals is not absurd; which
procedures accept a set, and how, is solved when the case arises.
`CTL_PLAN.md` section 11 is where that question continues.

## 9. Plan of attack

1. `loop/`: questions, kinds, reclassification, knowledge, the epoch as
   a rendez-vous with a rule policy; `prepare` rewritten as steps over
   them, behaviour identical, the three callers unchanged. The regression
   is the oracle checks of `Petri/test/`.
2. `proc/` with the adapters for what exists; `exam/` with the StateSpace
   and reachability rows; the walk driver becomes the loop of section 5 for
   those two families. Measured against the campaign baseline of
   2026-09-11.
3. `symbolic/`: hsc-pn's Petri net side moved here, libHSC as a linked
   library, `-hsc`; hsc-pn in libHSC becomes a thin example client or goes.
4. The total examinations and UpperBounds as rows, the conjoined goal and
   the isolation of subgoals measured on the QuasiLiveness flood; then
   CTL, with the checker as the dear procedure of the CTL kind.
5. The unfolder ported into `parse/` or a `coloured/` folder, then the
   skeleton decisions as procedures; the coloured strategies as rows.
6. ITS-Tools' `PetriSpotWalker` call sites replaced by one call per
   examination; the Java flows retired one by one; LTL last, on the
   semigroup engine.

Risks, named: the coloured unfolder has subtle semantics and the port must
be checked against the Java on the corpus; the SMT reductions' parity
(implicit places) needs the LP refiners and an exact checker before a rule
removes a place on a floating-point verdict, and a real solver behind
`lp/` when the refiners ask for one; and the loop's tuning is a campaign
question, not a design one, so every step above ends with a cluster run
read against the pages, honestly.
