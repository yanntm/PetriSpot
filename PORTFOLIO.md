# A portfolio that can be reasoned about

A design note for an overhaul of the solving loop of ITS-Tools + PetriSpot.
It says what the loop does today, what it measurably costs us, which objects
are missing, and — at the end — an ordered list of goals with the evidence each
one is supposed to move. Nothing here is implemented.

Everything numeric comes from the MCC 2026 campaign in
`Petri/test/mcc/csv/2026-09-05/`: 1953 instances, six examinations, 1800 s and
4 cores per test. The probes that produced the numbers are in
`Petri/test/probes/`.

## Why

The loop of structural reductions, SMT and random walks was written as a
front-end to the decision diagrams. It is now the engine: the walk answers 20 %
of the QuasiLiveness atoms, 59 % of the StableMarking ones, 42 % of the
UpperBounds ones, and the diagrams pick up what is left. Our results are state
of the art; the question is where the next verdict comes from.

Three measurements frame the answer.

* The walker received **90 h of a 1812 core-hour budget: 5 %**.
* In the 601 runs that burned the whole 1800 s and reported nothing, the run
  spends **48 % of its time on Liveness, 41 % on QuasiLiveness, 31 % on
  StableMarking** past the last timestamped line the Java side prints — a tail
  in which only the single threaded diagram engine is working
  (`probes/dd_tail.py`).
* The per property reduction phase of UpperBounds spent **47 319 s to produce
  4 502 s worth of verdicts**: 250 attempts, 68 of which closed something
  (`probes/perprop_payoff.py`).

None of those is a bug. They are the shape of a portfolio that cannot measure
itself and therefore cannot choose.

## What the loop already does

`ReachabilitySolver.applyReductions` is a fixpoint over a shrinking support,
not a pipeline. Per iteration:

1. build a `StructuralReduction` from the current net and a `RandomExplorer`
   over it, collect the still unsolved properties into `tocheck`;
2. walk — one `PetriSpotWalker.runReachability` request for the whole
   `tocheck` set, budget `30 + 5·min(|tocheck|, 50)` seconds, or the Java
   explorer's random / best-first / probabilistic cascade when PetriSpot is
   off;
3. SMT over the state equation, `smttime` escalating 5, 45,
   `45 + 15·iterations`, then each Parikh solution replayed as a guided walk;
4. `spn.computeSupport()`, `sr.setProtected(support)`, reduce again — with the
   SMT-based rules if the plain ones found nothing — and read the smaller net
   back;
5. `checkInInitial`, constant places, logic simplification, and on a barren
   iteration the atomic reducers.

A `verifyWithSDD` thread is started before the loop and runs beside all of it
on half the budget. So **the observation set is already first class**, under
the name *support*, recomputed every iteration so that each closed property
buys the reducer freedom; and **the diagrams already run concurrently** rather
than after.

Examinations reach that loop by decomposition.
`GlobalPropertySolver.buildQuasiLivenessProperty` emits one `EF ENABLED(t)` per
non-dominated transition — `computeDominatedTransitions` prunes first — after a
`ReductionType.LIVENESS` pass; StableMarking and OneSafe decompose the same
way; `applyExhaustiveMethods` hands the cohort to the loop with `timeout = -1`.

Per property models exist too, in `Application.java`, wherever the cohort is
the 16 formulas of an MCC examination: the CTL and LTL paths reduce each
property alone before solving, and UpperBounds does a global pass and then one
`UpperBoundsSolver.applyReductions` per property that is still open.

On coloured models the cheap abstractions run before the expensive ones on
purpose: the skeleton costs almost nothing and unfolding can exhaust memory, so
whatever the skeleton can settle is settled first. That ordering is right and
should survive the overhaul.

## What it does not do

**It ends on no progress, not on budget.** The condition is
`(iterations <= 1 || iter > 0) && properties remain`: one iteration that closes
nothing ends the loop, which drops through to the read-arc over-approximation
and the diagrams. Whether 100 s or 1500 s remain changes nothing. That, not an
absent walker, is why the deadlock runs gave the walk 13.7 s of 1800: it was
called, found nothing twice, and the loop concluded rather than reinvested.

**Effort is constants, not a resource.** `steps` is 10 000 then 1 000 000;
`smttime` is 5, 45, `45 + 15·iterations`; the walker budget is
`30 + 5·min(|tocheck|, 50)`. None is derived from the time left, from how many
cores are idle, or from what has been paying off on this instance.

**Nothing decides whether per property is worth it.** Where it exists it is
unconditional. On UpperBounds, 121 of 1953 logs reach the phase; those 121 run
250 per property reductions of which 68 close something — a 27 % hit rate for
**47 319 s of core time, 43 000 of them spent on nothing**. And when it works
it works narrowly: all 68 closed the very property they targeted, never a
neighbour. The wins are concentrated in a handful of families —
`HouseConstruction`, `DNAwalker`, `Murphy` — which is exactly the signal a
policy could use and no policy reads.

**Per property does not scale to the decomposed examinations.** QuasiLiveness
on a large instance is tens of thousands of atoms; one model per atom is not a
budget question, it is impossible. Any answer has to be *selective*, which
means it has to be *informed*.

**Support is a union, never a partition.** `sr.setProtected(support)` protects
what any remaining property observes, so two hard properties with disjoint
supports keep each other's places alive.

**Exploration state cannot survive an iteration.** Each pass builds a fresh
`RandomExplorer` over a freshly reduced net and spawns a fresh `petri64`;
UpperBounds alone made 10 087 such invocations, each re-serializing the net
through KERS and restarting from m₀, 301 of them killed by their timer with no
partial result. The loop is *right* to discard a frontier when the net changed
underneath it — but only because nothing maps the old markings onto the new
net.

**`DoneProperties` is a verdict map, not a knowledge base.** It records the
technique that closed a property. A bound that moved but did not close, a best
marking, a witness, a lease spent and lost: none of it is shared between
strategies or iterations, and none of it is recorded for us to tune against.
That asymmetry is why learning that the walker holds 5 % of the cores required
parsing 419 MB of logs.

## The objects

### Fact

| kind | content | who produces it |
| --- | --- | --- |
| verdict | formula φ is TRUE / FALSE | anyone |
| witness | a firing sequence reaching a marking satisfying φ | walks, BMC, diagrams |
| bound | expression *e* ∈ [lo, hi] over reachable markings | walks (lo), invariants and SMT (hi) |
| invariant | a P-flow or P-semiflow | PetriSpot |
| structure | *t* dead, *t* quasi-live, *p* constant, *p* implicit, a siphon, a trap | reductions, SMT, walks |
| exploration | best marking seen, states seen, depth reached | walks |
| profile | see below | everyone, as a side effect |

A witness is the most valuable kind: checkable, and it transfers. Making
PetriSpot excellent at producing witnesses is the same work as making the
portfolio's knowledge cheap to trust. Exploration facts are the least valuable
and the only ones that do not survive a change of model.

### Goal

A question with a knowledge state that strategies **narrow**: `unknown` /
`TRUE` / `FALSE` for a boolean, an interval `[lo, hi]` for a bound, closed when
`lo = hi`. The interval is not a metaphor — UpperBounds already runs this way,
with `Max Seen` rising from walks and `Max Struct` falling from invariants and
SMT, printed and forgotten between calls. On `RERS17pb113-PT-7` the log shows
`Max Seen` creeping from `[4,1,2,…]` to `[4,2,2,…]` while `Max Struct` sits at
7 for every expression: 92 s of walking that narrowed nothing, invisible to
anyone who might have changed strategy.

A goal also carries its **observation set**, its **provenance** (the goal it
was spawned from, and by which rule), and what has been spent against it.

### Model version

A net plus a **morphism back to its parent**. Versions form a DAG rooted at the
parsed net; an edge is a reduction, an unfolding, a skeleton abstraction or a
projection onto a smaller observation set. Each edge carries two maps or the
version is useless:

* **forward, on observations** — a parent place is a linear expression over the
  child's places, or a constant. This restates a bound or a formula on the
  child.
* **backward, on witnesses** — a child transition lifts to a sequence of parent
  transitions. This is how a trace found on a reduced net becomes a trace of
  the model the user asked about. The agglomeration rules already know that
  sequence and currently discard it.

### Profile

The cheap facts a run learns about the *model* rather than about a goal, and
the thing the current code has no place to put:

| signal | how it is obtained | what it should gate |
| --- | --- | --- |
| reduction yield | places removed by a pass / places before | whether to re-reduce, whether per goal models can pay |
| dead transition density | SMT dead-transition rules, walk coverage | whether the SMT reduction rules are worth their budget |
| walk productivity | verdicts per walk second, over iterations | whether to reinvest in walking or concede |
| state equation behaviour | UNSAT rate, solve time per problem | whether SMT deserves an escalating budget here |
| state space scale | invariant bounds, skeleton size, token counts | whether a diagram attempt is plausible at all |
| coloured structure | skeleton size versus unfolding estimate | whether to unfold, and when |

Reduction yield alone is a weak predictor and should not be used alone: over
3400 logs, runs that answered everything removed 38.9 % of places on the first
pass against 22.0 % for runs that left something unanswered, and 27.0 % versus
35.2 % removed nothing at all — a real difference, but the 0–1 % yield bucket
averages 0.42 unanswered values per log while the 1–10 % bucket averages 1.29
(`probes/reduction_resistance.py`). A profile has to be several signals, and
the framework's job is to make recording them free rather than to guess which
one matters.

### Strategy, cohort

A **strategy** is a worker with an applicability predicate and a cost profile:
random / best-first / Parikh-guided walk, structural reduction, invariants,
SMT, k-induction, BMC, decision diagrams, LTSmin. A **cohort** is the open
goals for one examination on one instance, the model DAG they live on, and the
budget they share. The cohort is the unit policy is set on: *this instance has
1800 s and 4 cores, spend them*.

## Observation sets, past the union

A reduction is legal exactly when it does not disturb what is observed, so the
observation set of the *open* goals bounds the reducer, and every goal that
closes makes the model smaller. The loop already exploits this and the
framework must keep it, not reinvent it. Two things the union cannot express.

**Per goal models, when earned.** A single goal observes almost nothing and its
private reduction can be drastic — that is why the phase exists and why it
pays 27 % of the time. Making it selective needs a predicate, and the profile
is where the predicate's inputs come from: yield on the shared model, whether
this goal's support is disjoint from the others', whether the family has paid
before. The measurement to beat is on record: 47 319 s spent, 4 502 s
productive.

**Priority inside the cohort.** The union treats every open goal as equally
worth protecting; a goal that has resisted five iterations costs the reducer
exactly what a fresh one does. There is no way to say "park this one, reduce
for the others, come back to it with a smaller model".

## What survives a reduction

The place a framework earns its keep. Each rule owes an answer to: which facts
about the parent are still facts about the child, and which facts about the
child lift to the parent?

* **verdicts on preserved formulas** lift child → parent; that is the rule's
  soundness argument, to be stated per rule rather than assumed;
* **witnesses** lift child → parent through the backward map;
* **invariants** project parent → child through the forward map when the rule
  is a place substitution, and are lost when it is not;
* **bounds** project when the observed expression survives the forward map;
* **structural facts** transfer for surviving objects, and a fact about a
  deleted object is not lost but promoted: "t is dead" is *why* the rule could
  delete t;
* **exploration facts do not survive** — a visited set, a frontier, a
  probability table: the state space changed underneath them.

That last line is what costs us. It argues for walking on the base model from
t = 0 rather than waiting, and for **projecting best markings** through the
forward map to reseed rather than restarting at m₀. Whether the projection of a
reachable marking is reachable in the child is a proof obligation per rule;
where it fails, the projected marking is still a useful heuristic target.

## The scheduler

A pool of (goal, model version) pairs. Workers own cores and **steal**: a free
worker takes the best item it can apply to, leases it, reports facts. Value is
a policy; a first version can be static, and the accounting will say what to
change. What must hold whatever the policy:

* **Leases are core-seconds**, not wall seconds. The only reason the 5 % figure
  was visible at all is that the walker prints its own duration.
* **A lease is interruptible**, and yields partial facts. Killing the walker's
  process is an interruption with no partial result: everything it learned dies
  with it, 301 times in this campaign.
* **A worker never blocks the pool.** The loop already overlaps one diagram
  thread with everything else; the stages *inside* an iteration serialize. Z3
  is barely multithreaded and PetriSpot's semiflow computation hit its 60 s
  wall 138 times, 101 of them on the `PolyORB*` family. Those are the moments
  the idle cores should be walking.
* **Facts are published, not returned.** One QuasiLiveness walk answers every
  transition it fires; publishing is what makes that free.

## Two mechanisms that need no framework

Both are small, both are testable against this campaign, and both can be built
before any of the above.

**The loop goes dormant instead of exiting.** When an iteration closes nothing,
leave a PetriSpot walking on the current model with the still-open goals and
let the run proceed to its next phase. If the walker reports a verdict, that is
new knowledge on the shared support, so retrigger a pass: a closed goal shrinks
the support, which unblocks reductions, which may unblock SMT. The loop's
termination becomes "no progress *and* the walker has nothing", instead of "no
progress this iteration".

**A PetriSpot alongside every diagram attempt.** When the run does leave the
loop and hands the model to the diagrams, spawn a PetriSpot on that same model
and goal, with a **memory budget** as well as a time budget, and let it run for
as long as the diagram engine does. It is free: the tail measured above is
348–867 s of single threaded diagram work per failed run, 48 % of a failed
Liveness run and 41 % of a failed QuasiLiveness one.

It should not be the same walker. A model that reached this point has resisted
everything cheap, so the strategies worth running are the ones that assume a
very large state space and know that a short random sweep will not do it:
restarts with diverse seeds, longer horizons, memory spent on a visited set to
avoid re-treading, directed search on the goal's support, exhaustive search of
a projection. That is a research direction of its own, and it is where "make
PetriSpot really good at finding traces" cashes out.

## A walk through QuasiLiveness

Today: dominated transitions pruned, a LIVENESS reduction pass, one
`EF ENABLED(t)` goal per surviving transition, then the reachability loop
iterating walk / SMT / reduce until an iteration closes nothing. What changes:

**t = 0.** Parse. On a coloured model the skeleton runs first, as now, and
closes what it can — a transition dead in the skeleton is dead in the
unfolding — because it is cheap and unfolding can exhaust memory.

**t = 0, same instant.** Walkers start on the unreduced net rather than after
the LIVENESS pass. Two cores, several threads, different strategies and seeds.
Every transition fired closes a goal; the knowledge base is a bitmap and
publishing is setting a bit. The gain is not the walk — the loop walks too — it
is that the walk overlaps the reduction instead of alternating with it, and
that its 30 s budget stops being a constant.

**t = ε, continuously.** The reducer keeps running against the current support
as it does now, but no longer waits for a walk to finish to see the support
move, and the walkers no longer restart at m₀ on each version: best markings
are projected forward.

**t = later.** The residue is the transitions no walk has fired. Today the loop
stops here. Instead, consult the profile: if reduction yield on this model is
good and the residue is small, spawn per goal models and hand them to SMT,
k-induction and the diagrams; if the model resists reduction — the case that
correlates with hardness — do not pay 47 000 seconds to learn it again, and
buy walk diversity with the remaining budget instead.

**Close.** QuasiLiveness is TRUE iff every child is TRUE and FALSE as soon as
one is FALSE, so the parent still closes early on the negative side. What is
new is that the answer carries what each strategy cost, not only which ones
contributed.

## What PetriSpot has to become

Today the interface is one shot: spawn `petri64`, hand it a KERS file and a
goal list, read stdout, kill it when its computed budget expires. That makes
four things impossible — retargeting mid-flight, reporting a partial bound,
keeping walk state across calls, being told the model just got smaller — and it
pays a full net serialization per call, 10 087 times in one examination.

The two shapes discussed are not alternatives:

* **A command language is the mechanism.** A long-lived `petri64` reading
  commands and streaming events: `load` a model, `derive` a version with its
  forward map, `goal add` / `goal drop` with priorities and budgets, `seed` a
  walk from a marking, `interrupt`; events `found <goal> <verdict> <witness>`,
  `bound <goal> <lo>`, `progress <goal> <stats>`. The libHSC surface syntax is
  the natural carrier — s-expressions, a vocabulary the two projects already
  share and what `Petri/test/logs/*.sexpr` already speaks.
* **A Java driver is the policy.** It holds the cohort, the model DAG, the
  profile and the budget, and decides when to interrupt. It must not hold the
  walk state.

Inside `petri64`, the same discipline as the outer scheduler: hot state thread
local, shared knowledge through one explicit documented structure. Two cores of
diverse walks means a small work-stealing pool over the goal set, with the
knowledge base — closed goals, best markings, the visited digest when memory is
granted — as the only shared object.

## Goals, in order

Each step names what it changes and the number it is meant to move. The order
is chosen so that every step is measurable on a rerun of this campaign before
the next one starts.

**G1 — Account for effort.** Every strategy reports core-seconds against a goal
and a model version, into a structured record the harness can read. No
behaviour changes.
*Moves*: makes every later step measurable; retires `Petri/test/mcc/`'s
log-scraping reconstruction of the same thing.

**G2 — Do not exit the loop with time on the clock.** Replace the no-progress
exit with a dormant walker that can retrigger a pass. Then co-schedule a
PetriSpot with every diagram attempt, with a memory budget.
*Moves*: the 415 timeouts and the 544 missed values; the 48 % / 41 % / 31 %
single-threaded tails.
*Risk*: retriggering can thrash; the dormant walker must hold a lease like
anyone else.

**G3 — A knowledge base instead of a verdict map.** `DoneProperties` grows into
facts: bounds with both ends, best markings, witnesses, and the model version
each fact is about. Strategies publish rather than return.
*Moves*: nothing on its own — it is the prerequisite for G4 and G5, and it is
what lets a bound survive an iteration.

**G4 — A long-lived PetriSpot.** The command language and interruption, driven
by the current loop, no scheduler yet. Kills the per-call serialization and the
per-call amnesia; lets a killed walk report what it had.
*Moves*: the 10 087 invocations and the 301 killed calls; first place where
"good at finding traces" can be worked on directly.

**G5 — A model profile, and per goal models on evidence.** Record the profile
signals; gate the per property phase on them instead of running it always.
*Moves*: the 43 000 wasted seconds of the UpperBounds per property phase,
without giving up the 68 verdicts it does find.

**G6 — Versioned models with their maps.** The forward map on observations and
the backward map on witnesses, per reduction rule, with the survival table.
*Moves*: lets walk state be projected instead of discarded, and lets a witness
found on a reduced net be reported on the original.

**G7 — The scheduler.** Work stealing over (goal, model version), leases,
priorities. Last, because only G1 gives the evidence to argue about its policy.

## Open questions

1. Which facts survive which reduction rule, in which direction? A table, rule
   by rule. Prerequisite for G6, not a refinement of it.
2. Is the projection of a reachable marking reachable in the reduced net? Per
   rule. Where not, is it still a good target?
3. What does lifting a witness through a long chain of versions cost, and do we
   keep the chain or compose eagerly?
4. How many model versions can live at once before memory decides for us?
5. What is the value function for stealing, and can it be learned from the
   campaign record rather than guessed?
6. When the cohort is 33 676 goals, every per-goal object must be an index into
   a flat array. What is the smallest goal representation that still carries
   provenance?
7. What does a walker do with a memory budget that it cannot do without one,
   and how much does a visited digest actually buy on the models that reach the
   diagram phase?
8. Where does the roll-up rule live — the one turning child verdicts into the
   parent's — and how much is shared between QuasiLiveness, StableMarking and
   OneSafe?

The numbers to beat are in `Petri/test/mcc/csv/2026-09-05/`: 544 values the
consensus has and we do not, 1237 nobody has, 415 runs that spent 1800 s to say
nothing, and a walker holding 5 % of the cores.
