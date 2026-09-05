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
* Rerunning ReachabilityDeadlock across the arrival of `walk/DeadlockStrategy.h`
  gained 19 instances and 8 values the consensus lacks, and lost 8 — and the
  eight are mostly not lost to a hard model: rerun twice, one turns out to be an
  infrastructure flake and four are found intermittently within 65 s of an
  1800 s budget, always on the second walker call, while the failing runs give
  up after about a second and spend 1796 s in `its-ctl`
  (`csv/2026-09-06-RD/`).

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
UpperBounds alone made 10 087 such invocations, each re-serializing the net and
restarting from m₀. The verdicts survive — `PetriSpotWalker` parses stdout as
it arrives, so the 301 calls killed on their timer kept what they had already
printed — but the exploration state has no output line and dies every time. The
loop is *right* to discard a frontier when the net changed underneath it; the
waste is that it is discarded even when the net did not.

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

### Model, and what a fact is about

A full DAG of model versions with maps both ways is a large piece of machinery
and it is **not** on the critical path; it is parked here so the rest can be
written without it.

What is needed instead is the little that already works, made explicit. Place
names survive the transformations — that is how the code reindexes after a
structural reduction — so a fact stated over places transfers by name.
Transition names do not: on large nets they are mangled rather than composed,
because composing them blows up in size, which is precisely why a *witness*
cannot be lifted today and why lifting is parked with the DAG.

So a fact is only useful if it names what it is about, and the two things worth
naming now are:

* **the base model**, for anything proven of the net itself — a bound, an
  invariant, a verdict on a property of the original question;
* **the generation**, a counter bumped every time the reduction loop produces a
  new net, for anything that is only meaningful on the net that produced it —
  a best marking, a firing profile, an exploration statistic. A fact from an
  older generation is not wrong, it is stale, and the cheap policy is to drop
  it rather than to translate it.

One tag more is worth carrying, because the loop already builds such models and
nothing records which kind they are: a derived net is an **over-approximation**
(the read-arc abstraction, the skeleton), an **under-approximation**, or
**exact**. A verdict is conclusive in one direction only for the first two, and
today that direction is re-derived at each use site.

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

Parked with the model DAG, and recorded here because it is the obligation any
later attempt inherits: each rule owes an answer to which facts about the
parent are still facts about the child, and which facts about the child lift to
the parent. Verdicts on preserved formulas lift child to parent — that is the
rule's soundness argument. Invariants and bounds project the other way when the
rule is a place substitution. A structural fact about a deleted object is not
lost but promoted: "t is dead" is *why* the rule could delete it.

The one line that matters before any of that is built: **exploration facts do
not survive** — a visited set, a frontier, a firing profile. The state space
changed underneath them. Today that costs nothing extra, because exploration
state does not survive the *call* either; that is what makes a generation
counter enough for now, and what makes G4 worth doing before G6.

## A pool, not a market

"Work stealing" was too strong. What is wanted is smaller and buildable: **a
pool of work for a walker, a budget attached to each item, and control over
what we handed out.** The portfolio decides; the walker executes what it is
given and reports what it spent.

That means, concretely:

* **an item is (goals, model, budget)** — the goals a walker should try, the
  net they are stated on, and what it may spend. Today the budget is
  `30 + 5·min(|tocheck|, 50)` seconds and nothing else is passed;
* **the caller can revise an item** — drop goals that closed elsewhere, raise
  the budget of one that looks promising, cancel one whose model went stale;
* **effort is reported back** per item, which is what makes the profile
  possible at all.

Triggers already exist in the loop and are the model for the rest: if a pass
found dead transitions, try again; strange arc values trigger the SMT dead
transition test. They are hand-written conditions on facts, which is the right
shape — the missing part is a place to keep the facts they read, and a record
of what each trigger cost when it fired.

Two rules to hold whatever the policy:

* **Effort is core-seconds**, reported per item, not wall time. The only reason
  the 5 % figure was visible at all is that the walker prints its own duration.
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

Both are instances of one rule: **a strategy that finds nothing quickly should
reinvest, not conclude.** The deadlock rerun shows the same failure inside
PetriSpot, where `WalkDriver.h` ends its rounds after about a second of a
thirty-five second budget, so the rule applies at both levels.

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

Two claims I made for a long-lived interactive `petri64` do not survive
reading the code, and the case has to be made without them.

`PetriSpotWalker.run` starts a thread that parses stdout *as it arrives*, so a
walker killed on its timer keeps every verdict and every bound it had already
printed — the counts are taken after the kill. **Partial results already
survive.** And the spawn is not the expense either: the walker's own reports
are hundreds of milliseconds, of which the net serialization is a small part.
The 10 087 invocations are ugly, not costly.

What actually dies with the process is the thing that has no output line: the
**exploration state**. Firing counts accumulated over every walk, which regions
of the net turned out hard to enter, which restart states were promising, which
heuristic was fast here and which was slow. That is exactly the content that
would make a walker good at models with a very large state space, and it is
recomputed from nothing on every call.

So the right shape is not an interactive protocol — which is race-prone,
error-prone and buys retargeting we have no policy to use yet — but **fewer,
bigger calls**:

* one call per (net, goal set, budget), given the *whole* budget the portfolio
  is willing to spend now, not a sliver of it;
* all the diversity and all the memory **inside** the process: restarts drawing
  from a shared pool, firing profiles kept across restarts, rare-event
  targeting, self-evaluation of which strategy is paying, focus on the parts of
  the system that are hard to enter;
* results streamed as they are found, as now, so a kill costs only the memory
  and not the verdicts.

The evidence that this is where the money is sits in the deadlock rerun. On
`ShieldPPPt-PT-040B` the walker is called six times and returns in 1045, 643,
337, 1146, 409 and 371 ms against a 35 s budget, each time having solved
nothing, because `WalkDriver.h` stops the rounds when one "solved nothing,
improved no bound, and every walk in it ended on the step budget" — about a
second, on a single property. The run then spends 1796 s in `its-ctl`. Thirty
four of thirty five seconds unspent, then half an hour spent elsewhere.

A command language remains the right *mechanism* the day a policy needs it, and
the libHSC s-expression surface is its natural carrier. It is not the next
step.

## Goals, in order

**G1 — A record of effort, in `DoneProperties`.** Not more tracing: the same
structure, extended. It already is the portfolio's synchronisation point, so it
grows the *open* side of the ledger — the properties still to prove and the net
generation they are being worked on — plus, per task, the time it took and a
measure of what it returned for that time. Default cost is a few counters; a
verbose mode may justify every rule applied, and is off by default because at
scale that is thousands of models and traces.
*Moves*: nothing by itself. It is what makes G2 and G5 arguable instead of
guessed, and it retires `Petri/test/mcc/`'s reconstruction of the same numbers
from 419 MB of logs.

**G2 — Budgets from the clock, not from constants.** `steps`, `smttime` and the
walker budget become functions of the time left and of what the last pass
returned. Same rule one level down, inside PetriSpot: the round scheduler
should not concede after a second of a thirty-five second budget.
*Moves*: the eight deadlock instances lost to the round rule; the 415 timeouts;
the 48 % / 41 % / 31 % single-threaded tails.
*Risk*: spending more on a strategy that will not find anything. G1 first, so
the escalation can be read afterwards.

**G3 — Listeners on that record.** Strategies subscribe to it: a closed goal
interrupts the walker still working on it, a new dead transition retriggers the
reduction pass, a bound that moved reopens a comparison. The triggers already
in the loop — try again if a pass found dead transitions, run the SMT dead
transition test on strange arc values — become subscriptions rather than
open-coded conditions.
*Bounded on purpose*: facts, counters and a generation number. No model DAG, no
witness lifting, no per-rule justification unless asked for.

**G4 — Fewer, bigger PetriSpot calls, with memory inside.** One call per (net,
goals, budget) carrying the whole budget; firing profiles, restart pools,
rare-event targeting and strategy self-evaluation kept across restarts within
the call; results streamed as now.
*Moves*: the models that need a large-state-space search rather than a sweep —
where "really good at finding traces" is actually won.

**G5 — Per goal models on evidence.** Gate the per property phase on the
profile accumulated in G1 and G3 instead of running it unconditionally.
*Moves*: the 43 000 wasted seconds of the UpperBounds per property phase,
without giving up the 68 verdicts it does find.

**G6 — Versioned models with maps both ways.** Deprioritised. Needed the day a
witness must be reported on the original net, or walk state must cross a
reduction. Until then the generation counter and the place-name identity carry
what we use.

**G7 — Scheduling policy.** Start small: a pool, budgets, and the ability to
revise an item. Revisit once G1 says where the time goes.

## Decided

* **A fact must be usable, not merely true.** Anything proven of the net is
  stated over the base model; place names carry across transformations, which
  is how reindexing already works. Transition names are mangled on large nets,
  so a firing sequence is not liftable today — the reason witness lifting waits
  for G6.
* **A derived net carries its direction**: over-approximation, under-approximation
  or exact. The loop builds all three and re-derives the direction at each use.
* **Not many versions at once**, so the memory question does not arise; a
  generation counter is enough to tell stale exploration state from fresh.
* **A pool, not a market.** Items are (goals, net, budget), the caller keeps
  the right to revise them, and the walker reports what it spent.
* **Below a threshold of open goals the regime changes.** With thousands of
  goals, per-goal bookkeeping and per-goal models are both out; the framework
  has to work in bulk first and switch to individual attention when the residue
  is small. Where that threshold sits is a measurement, not a constant to pick
  now.
* **The three global examinations share one essence**: a large set of queries
  to eat through, with an early break — one FALSE settles QuasiLiveness, one
  unstable place settles StableMarking, one unsafe place settles OneSafe. The
  roll-up rule is the same object three times.
* **A walker with memory** is one that keeps firing counts over all its walks,
  aims at rare events, studies the shape of the state space it is in, reuses
  restart states, remembers which heuristics ran fast or slow here, and focuses
  on the parts of the system that are hard to get into.

## Still open

1. Where the "few enough goals" threshold sits, per examination.
2. What the efficiency measure of a task should be, given that most tasks
   legitimately return nothing. Verdicts per core-second is too coarse; bound
   movement and support shrinkage are candidates.
3. Whether the round rule in `WalkDriver.h` should be relaxed, replaced by an
   escalation, or made conditional on the budget left. Rerunning the eight lost
   deadlock instances twice says four of them are found intermittently, always
   on the second walker call and within 65 s of an 1800 s budget, so a third
   call is the experiment. Three are lost consistently and are a heuristic
   question, not a budget one.
4. How much a visited digest actually buys on the models that reach the diagram
   phase, and at what memory.

The numbers to beat are in `Petri/test/mcc/csv/2026-09-05/` and
`csv/2026-09-06-RD/`: 544 values the consensus has and we do not, 1237 nobody
has, 415 runs that spent 1800 s to say nothing, a walker holding 5 % of the
cores, and eight deadlock instances lost to a one second give-up.
