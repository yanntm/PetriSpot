# A portfolio that can be reasoned about

This is a design note, not a description of what exists. It proposes a shape
for the solving loop of ITS-Tools + PetriSpot: what a goal is, what we know at
any moment, which model that knowledge is about, and who gets a core. Nothing
here is implemented. The point is to have something to argue with.

## Why

The loop of structural reductions, SMT and random walks was written as a
front-end to the decision diagrams. It is now the engine: on the MCC 2026
campaign of `Petri/test/mcc/csv/2026-09-05/`, the walk answers 20 % of the
QuasiLiveness atoms, 59 % of the StableMarking ones, 42 % of the UpperBounds
ones, and the diagrams pick up what is left. The intent of the architecture and
its actual behaviour have drifted apart, and the accounting shows where that
hurts.

The walker received **90 h of a 1812 core-hour budget: 5 %**. In the 601 runs
that burned the whole 1800 s and reported nothing, the mean time handed to the
walker was 13.7 s for ReachabilityDeadlock, 116 s for Liveness, 282 s for
UpperBounds. Nobody decided that. The next section says what the loop actually
does, because the shape of the answer depends on it.

## What the loop already does

`ReachabilitySolver.applyReductions` is the loop, and it is a fixpoint over a
shrinking support, not a pipeline. Per iteration:

1. build a `StructuralReduction` from the current net and a `RandomExplorer`
   over it, collect the still unsolved properties into `tocheck`;
2. walk — one `PetriSpotWalker.runReachability` request for the whole
   `tocheck` set, budget `30 + 5·min(|tocheck|, 50)` seconds, or the Java
   explorer's random / best-first / probabilistic cascade when PetriSpot is
   off;
3. SMT over the state equation, with `smttime` escalating 5, 45,
   `45 + 15·iterations`, then replay each Parikh solution as a guided walk;
4. `spn.computeSupport()`, `sr.setProtected(support)`, reduce again — with
   SMT-based rules if the plain ones found nothing — and read the smaller net
   back;
5. `checkInInitial`, constant places, logic simplification, and on a barren
   iteration the atomic reducers.

A `verifyWithSDD` thread is started before the loop and runs beside all of it
with half the budget. So the two things a first reading of the logs suggests
are missing are in fact there: **the observation set is already first class**,
under the name *support*, recomputed every iteration so that each closed
property gives the reducer more freedom; and **the diagrams already run
concurrently** with the whole loop rather than after it.

The examinations reach that loop by decomposition. `buildQuasiLivenessProperty`
emits one `EF ENABLED(t)` per non-dominated transition — `computeDominatedTransitions`
prunes first — after a `ReductionType.LIVENESS` pass, and `applyExhaustiveMethods`
hands the cohort to the loop with `timeout = -1`. StableMarking and OneSafe
decompose the same way.

## What it does not do

* **The loop ends on no progress, not on budget.** Its condition is
  `(iterations <= 1 || iter > 0) && properties remain`: an iteration that
  closes nothing ends the loop, which then falls through to the read-arc
  over-approximation and the diagrams. Whether 100 s or 1500 s of the budget
  are left changes nothing. That, not an absent walker, is why the walk got
  13.7 s of 1800 in the deadlock runs: it was called, it found nothing twice,
  and the loop concluded rather than reinvested.
* **Effort is a constant, not a resource.** `steps` is 10 000 then 1 000 000,
  `smttime` is 5 then 45 then `45 + 15·iterations`, the walker budget is
  `30 + 5·min(|tocheck|, 50)`. None of them is derived from what remains of the
  budget, from how many cores are idle, or from which strategy has been paying
  off on this instance.
* **Exploration state cannot survive an iteration.** Each pass builds a fresh
  `RandomExplorer` over a freshly reduced net and spawns a fresh `petri64`;
  UpperBounds alone made 10 087 such invocations, each re-serializing the net
  through KERS and restarting from m₀, and 301 of them died on their kill
  timer with no partial result reported. The loop is *right* to discard a
  frontier when the net changes underneath it — but only because nothing maps
  the old markings onto the new net.
* **Support is a union, never a partition.** `sr.setProtected(support)` protects
  what *any* remaining property observes. Two hard properties with disjoint
  supports keep each other's places alive, and no per-property model is ever
  built.
* **`DoneProperties` is a verdict map, not a knowledge base.** It records a
  verdict and a technique string per property. A bound that moved but did not
  close, a best marking, a witness, a lease that was spent and failed: none of
  it is shared between strategies or between iterations, and none of it is
  recorded for us to tune against.

## The objects

### Fact

The unit of knowledge. A fact is about a **model version** and is either

| kind | content | who produces it |
| --- | --- | --- |
| verdict | formula φ is TRUE / FALSE | anyone |
| witness | a firing sequence reaching a marking satisfying φ | walks, BMC, diagrams |
| bound | expression *e* ∈ [lo, hi] over reachable markings | walks (lo), invariants and SMT (hi) |
| invariant | a P-flow or P-semiflow | PetriSpot |
| structure | *t* is dead, *t* is quasi-live, *p* is constant, *p* is implicit, a siphon, a trap | reductions, SMT, walks |
| exploration | best marking seen, states seen, depth reached | walks |

A witness is the most valuable kind, because it is *checkable* and it
*transfers*: replaying it costs nothing and convinces anyone. Making PetriSpot
excellent at producing witnesses is the same thing as making the portfolio's
knowledge cheap to trust. The exploration facts are the least valuable and the
only ones that do not survive a change of model.

### Goal

A question with a knowledge state that strategies **narrow**:

* boolean goal — state is `unknown` / `TRUE` / `FALSE`;
* numeric goal — state is an interval `[lo, hi]`, closed when `lo = hi`;
* structural goal — a set to be decided element by element, closed when empty.

The interval is not a metaphor: UpperBounds already runs this way, with
`Max Seen` moving up from walks and `Max Struct` moving down from invariants
and SMT, and the current code prints both and forgets them between calls. On
`RERS17pb113-PT-7` the log shows `Max Seen` creeping from `[4,1,2,…]` to
`[4,2,2,…]` while `Max Struct` sits at 7 for every expression: 92 s of walking
that narrowed nothing, invisible to whoever might have changed strategy.

A goal also carries:

* **its observation set** — the places and transitions its answer depends on.
  This is the single most important field, because it is what decides which
  reductions are legal (below);
* **provenance** — the goal it was spawned from and the rule that spawned it;
* **a budget** and what has been spent against it.

### Model version

A net, plus a **morphism back to its parent**. Model versions form a DAG rooted
at the parsed net. An edge is a reduction, an unfolding, a skeleton
abstraction, or a projection onto a smaller observation set.

Each edge must carry two maps, or the version is useless:

* **forward, on observations**: a place of the parent is a linear expression
  over the child's places, or a constant. This is how a bound or a formula is
  restated on the child.
* **backward, on witnesses**: a transition of the child lifts to a sequence of
  parent transitions. This is how a trace found on a reduced net becomes a
  trace of the model the user asked about. Agglomeration rules already know
  this sequence; they currently do not keep it.

### Strategy

A worker with an **applicability predicate** and a cost profile. Random walk,
best-first walk, Parikh guided walk, structural reduction, invariant
computation, SMT over the state equation, k-induction, BMC, decision diagrams,
LTSmin. Applicability is a cheap test on (goal, model version): SMT wants a
linear property, the diagrams want a net under some size, the skeleton wants a
coloured model, a walk wants an existential goal.

### Cohort

The set of open goals for one examination on one instance, with the model DAG
they live on and the budget they share. The cohort is the unit the user sets
policy on: *this instance has 1800 s and 4 cores, spend them*.

## Observation sets, past the union

A reduction is legal exactly when it does not disturb what is observed, so the
observation set of the *open* goals bounds the reducer, and every goal that
closes makes the model smaller. The loop already exploits this: the walk feeds
the reducer through `computeSupport`, and the reducer feeds the walk back a
smaller net. That is the good part of what exists and the framework must keep
it, not reinvent it.

Two things the union of supports cannot express.

**Per-goal models.** `setProtected` takes the union, so two hard goals with
disjoint supports protect each other's places forever. Taken to its limit the
same argument gives one model per goal: a single goal observes almost nothing
and its private reduction can be drastic. That is where the expensive
strategies belong. The cost is that the model DAG fans out, so it is worth
doing only when the residual cohort is small — a policy question with a number
in it, to be measured, not assumed.

**Priority inside the cohort.** The union treats every open goal as equally
worth protecting. A goal that has resisted five iterations and a goal spawned
one iteration ago cost the reducer the same, and the loop has no way to say
"park this one, reduce for the others, come back with a fresh model".

## What survives a reduction

This is the hard part and the place where a framework earns its keep. Each
reduction rule owes an answer to: which facts about the parent are still facts
about the child, and which facts about the child lift to the parent?

Sketch of the expected shape:

* **verdicts on preserved formulas** lift child → parent. That is the soundness
  argument of the rule itself; it must be stated per rule, not assumed.
* **witnesses** lift child → parent through the backward map, which is why the
  backward map must exist.
* **invariants** project parent → child through the forward map when the rule
  is a place substitution, and are lost when it is not.
* **bounds** project when the observed expression survives the forward map.
* **structural facts** transfer for the objects that survive, and a fact about
  an object the rule deleted is not lost but *promoted*: "t is dead" is why the
  rule could delete t.
* **exploration facts do not survive.** A visited set, a search frontier, a
  probability table: the state space changed underneath them.

That last line is the one that costs us. It argues for two things at once:
start walking on the base model at t = 0 rather than waiting for the reducer,
because time spent waiting is not recoverable; and, when a new model version
appears, **project the best markings** through the forward map and reseed the
walk from them instead of restarting at m₀. Whether the projection of a
reachable marking is reachable in the child is a proof obligation per rule, and
where it fails, the projected marking is still a useful heuristic target even
if not a legal start state.

## The scheduler

A pool of (goal, model version) pairs. Workers own cores and **steal**: a free
worker scans the pool for the item with the best expected value it can apply
to, takes a lease on it, and reports facts back. Value is a policy; a first
version can be static — witness-shaped goals to the walkers, arithmetic-shaped
goals to SMT, everything to the reducer's observation-set watcher — and the
accounting will say what to change.

What the framework must enforce, whatever the policy:

* **Leases are core-seconds.** Not wall seconds. The reason we could see the
  5 % figure at all is that the walker prints its own duration; nothing else
  in the pipeline does.
* **A lease is interruptible.** A strategy holding a lease must poll a flag and
  yield promptly with whatever partial facts it has. Killing the walker's
  process is an interruption with no partial result: everything it learned dies
  with it, 301 times in this campaign.
* **A worker never blocks the pool.** The loop already overlaps one decision
  diagram thread with everything else; the stages *inside* an iteration are
  what serialize. SMT and the invariant solver are the cases that matter: Z3 is
  barely multithreaded, and PetriSpot's semiflow computation hit its 60 s wall
  138 times in this campaign, 101 of them on the `PolyORB*` family. Those are
  exactly the moments when the idle cores should be walking.
* **Facts are published, not returned.** A strategy that finds a witness for
  one goal often answers others by accident — a single QuasiLiveness walk
  answers every transition it fires. Publishing into the knowledge base is what
  makes that free.

## Traceability

`DoneProperties` records the technique that closed a property, which is what
feeds the TECHNIQUES field, and nothing else: not what was tried and failed,
not what it cost, not on which model version. That asymmetry is why the only
way to learn that the walker holds 5 % of the cores was to parse 419 MB of
logs. For each goal the record should carry the chain of model versions it was
attacked on, the leases spent on each, the facts that closed it, and the spawn
edges from its parent. Rolled up per instance that is an honest TECHNIQUES
field; per campaign it is the table saying where the 1800 s went. The collector
in `Petri/test/mcc/` reconstructs a shadow of it from log lines; it should be
reading it from the tool.

## A walk through QuasiLiveness

The examination asks whether every transition can fire in some reachable state.
Today: dominated transitions are pruned, a LIVENESS reduction pass runs, one
`EF ENABLED(t)` goal per surviving transition goes into the cohort, and the
reachability loop iterates walk / SMT / reduce until an iteration closes
nothing. What would change, step by step:

**t = 0.** Parse. If the model is coloured the skeleton is available before the
unfolding is, and it is an abstraction: a transition dead in the skeleton is
dead in the unfolding, so the skeleton closes negative goals for free while the
unfolder is still running. That is concurrency the current sequence cannot
express, since it needs the SPN before it builds any property.

**t = 0, same instant.** Walkers start on the unreduced net rather than after
the LIVENESS pass. Two cores, several threads, different strategies and seeds.
Every transition fired closes a goal; the knowledge base is a bitmap and
publishing a fact is setting a bit. The gain is not the walk itself — the loop
walks too — it is that the walk overlaps the reduction instead of alternating
with it, and that its 30 s budget stops being a constant.

**t = ε, continuously.** The reducer keeps running against the current support,
as it does now, but it no longer has to wait for a walk to finish to see the
support move, and the walkers no longer restart at m₀ on each new version:
their best markings are projected forward through the version's map.

**t = later.** The residual cohort is the transitions no walk has fired. The
loop today stops here, on its no-progress condition. Instead: if the residue is
small, spawn per-goal models — each observing one transition, so far more
reducible than their union — and hand those to SMT, k-induction and the
diagrams. If it is large, the instance is hard and the remaining budget buys
more walk diversity rather than one deep exhaustive attempt.

**Close.** QuasiLiveness is TRUE iff every child is TRUE and FALSE as soon as
one is FALSE, so the parent closes early on the negative side, as it already
does. What is new is that the answer carries what each strategy cost, not only
which ones contributed.

## What PetriSpot has to become

Today the interface is one shot: spawn `petri64`, hand it a KERS file and a
goal list, read stdout, kill it when its computed budget runs out. That
interface makes four things impossible — retargeting mid-flight, reporting a partial bound, keeping walk
state across calls, and being told that the model just got smaller — and it
pays a full net serialization per call, 10 087 times in one examination.

Two shapes were suggested. They are not alternatives:

* **A command language** is the mechanism. A long-lived `petri64` reading
  commands and streaming events: `load` a model, `derive` a version with its
  forward map, `goal add` / `goal drop` with priorities, `seed` a walk from a
  marking, `interrupt`, and events `found <goal> <verdict> <witness>`,
  `bound <goal> <lo>`, `progress`. The surface syntax of libHSC is the natural
  carrier — s-expressions, a vocabulary the two projects already share, and
  already what `Petri/test/logs/*.sexpr` speaks.
* **A Java driver** is the policy. It holds the cohort, the model DAG and the
  budget, and it is the thing that decides to interrupt. It must not hold the
  walk state.

Inside `petri64`, the same discipline as the outer scheduler and the same rule
as the rest of the code: hot state thread-local, shared knowledge through one
explicit documented structure. Two cores of diverse walks means a small work
stealing pool over the goal set, with the knowledge base — closed goals, best
markings per goal — as the only shared object.

## Open questions

1. Which facts survive which reduction rule, in which direction? This wants a
   table, rule by rule, and it is a prerequisite, not a refinement.
2. Is the projection of a reachable marking reachable in the reduced net? Per
   rule. Where it is not, is the projected marking still a good target?
3. What does a witness cost to lift through a long chain of versions, and do we
   keep the chain or compose the maps eagerly?
4. How many model versions can live at once before memory decides for us? The
   sparse structures are shareable, the question is what fraction is shared.
5. What is the value function for stealing, and can it be learned from the
   campaign record rather than guessed?
6. When the cohort is 33 676 goals, every per-goal object must be an index into
   a flat array. What is the smallest goal representation that still carries
   provenance?
7. Where does the examination's own semantics live — the roll-up rule that
   turns child verdicts into the parent's — and how much of it is shared
   between QuasiLiveness, StableMarking and OneSafe?

## A first increment

Nothing above needs to be built at once, and the useful order is the one that
produces a measurement at each step:

1. **Accounting.** Make every strategy report core-seconds against a goal.
   Costs nothing, changes nothing, and turns every campaign into evidence.
2. **A long-lived PetriSpot** with the command language and interruption, still
   driven by the current pipeline. Kills the per-call serialization and the
   50 s amnesia; measurable on the UpperBounds calls alone.
3. **A budget instead of a no-progress exit.** The loop's termination
   condition is one line; making a barren iteration reinvest the remaining
   time rather than concede is the smallest change with a measurable effect,
   and the 601 runs that ended early with time to spare are the sample.
4. **Walk concurrently with the reducer**, on QuasiLiveness first, with best
   markings projected across versions. This is where the model DAG and its
   maps become unavoidable.
5. **The scheduler**, last, once there is enough accounting to argue about its
   policy.

The numbers to beat are in `Petri/test/mcc/csv/2026-09-05/`: 544 values the
consensus has and we do not, 1237 nobody has, 415 runs that spent 1800 s to say
nothing, and a walker holding 5 % of the cores.
