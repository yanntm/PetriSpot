# `walk/` — explicit walk: data structures and the step

Everything here is sparse and per thread except `WalkNet`, which is built once
from the net and only read afterwards.

## WalkNet (shared, read-only)

* `pre(t)`: the column of `flowPT`, `(place, weight)` pairs sorted by place.
* `effect(t)`: `post(t) - pre(t)`, sorted by place, zeros removed. A transition
  with an empty effect is a self-loop: firing it never changes the marking.
* `consumersOf(p)`: the transitions with an input arc from `p`, each with its
  arc weight, **sorted by weight**. This is the place-to-transition matrix,
  the only cross-index we build. No transition-to-transition structure is ever
  computed.
* `maxWeight(p)`: the largest input arc weight among consumers of `p`.

## Marking (per thread)

A `SparseArray<T>` of the places with a non-zero token count. `apply(effect,
onChange)` merges the effect in place: for each `(p, d)` of the effect, the
entry of `p` is located by binary search, updated (`old + d`), inserted when it
was absent, removed when it reaches zero, and `onChange(p, old, new)` is
called. Cost `O(|effect| * log |marking|)` plus the tail shift of an insert or
a removal. Copy assignment between markings reuses the destination capacity,
so a reset to the initial marking does not allocate.

## EnabledSet (per thread)

* `unsat[t]`: number of input arcs `(p, w)` of `t` with `m[p] < w`. `t` is
  enabled iff `unsat[t] == 0`.
* `list`: the enabled transitions, packed; `position[t]` is the index of `t`
  in `list` or NONE. Removal swaps with the last element; random choice is
  `list[rand]`.

`initialize(marking)` computes `unsat` from scratch once (`O(arcs)`); the
walker keeps that state as a snapshot and restores it by copy at every reset.

`onPlaceChanged(p, old, new)` is the delta step, called by `Marking::apply`:

* `new < old`: consumers of `p` with weight `w` in `(new, old]` lose one
  satisfied arc: `unsat[t]++`, and `t` leaves `list` when it reaches 1.
* `new > old`: consumers with `w` in `(old, new]` gain one: `unsat[t]--`,
  and `t` enters `list` when it reaches 0.

Since `consumersOf(p)` is sorted by weight the affected range is found by two
binary searches and only the transitions that actually flip are visited. If
`old` and `new` are both at least `maxWeight(p)`, or both below the smallest
weight, nothing is visited. Cost per fired transition:
`O(sum over touched p of log(fanout(p)) + number of flips)`, independent of
the number of transitions or places.

## Target and TargetSet

A `Target` is one goal: a normal-form `Expression` evaluated on the marking
(`reached` is `goal.eval(marking)`), deadlock (`reached` is "the enabled list
is empty"), or a bound: a weighted sum of places to maximise, with an optional
limit (a known upper bound) whose value, once reached, ends the target as a
goal would (`reached` is `value >= limit`).

A `TargetSet` is every open goal of a run, shared by the threads, so that one
walk answers many questions:

* per target: name, MCC verdict word, the `Target`, and an atomic `solved`
  flag; `claim(i)` is a compare-and-swap, so exactly one thread wins a target
  and the count of open targets is exact;
* `mentions(i)`: the places target `i` mentions, each with the direction of
  change that can make `i` hold: an increase for a place with a positive
  coefficient in a `>=` atom (or a negative one in a `<=`), a decrease for the
  reverse, either for `==`, `!=` or a place appearing both ways; a bound grows
  like a `>=`. The index built from it is per walker (below);
* `deadlocks()`: the deadlock targets, checked only when the enabled list is
  empty;
* per bound target, an atomic running maximum: a walker publishes a value
  only when it beats the last one it published (a thread-local array), so
  the atomic is touched at improvements only; the driver prints it at exit.

Each walker owns its index: `own`, the targets it checks, and per place the
lists `up` and `down` of the own targets an increase or a decrease of the
place can make hold. A walker re-evaluates a target only when the fired
transition changed a place the goal mentions in the direction that matters:
the places touched by `Marking::apply` are collected with the sign of their
change, the `up` lists of the places that grew and the `down` lists of those
that shrank are merged with an epoch stamp per target (no clearing), and each
candidate still open is evaluated in full (goal expressions are small trees).
A `fireable` target, all lower bounds, is thus never checked on the places a
firing drains. Solved ids met while scanning a list are swapped out of it, so
the lists shrink over the run without a compaction pass.

With `--partition=N` a sweep over at least `N` targets is split between the
threads, thread `k` owning the ids congruent to `k`: the checks per step
divide by the thread count, at the price of the finds a thread walks past
because another owns them. Off by default: on ResIsolation-PT-N10P4 the split
doubled the step rate and claimed 30 % fewer targets in 20 s. A walker whose
own targets are all solved takes over every open target (`refill`).
The cost per step is proportional to the arcs touched and to the targets that
depend on them, never to the number of targets. Every open target is evaluated
in full once at the start of a run and after a reset to a pooled state; a
reset to the initial marking checks nothing, the first check already covered
it and no target can have started to hold there since. The proportionality
breaks on hub places: a place consumed by tens of thousands of transitions
puts all their `fireable` targets in its list, and a step that moves it
re-evaluates them all (WALK_PLAN.md section 9).

A walker has an optional *focus*: the target its strategy is built for. It
claims any open target it happens to reach, and stops when its focus is
solved (by itself or by another thread), when no target is open, or when the
budget is exhausted. A walker without focus (random strategy) runs until no
target is open or the budget ends: that is the sweep round of the driver.

## Strategy

`choose(context) -> transition` picks the next transition among the enabled
ones, or returns `RESTART` to ask the walker for a fresh run. The context
exposes the marking, the enabled set, the net and the RNG.

* `RandomStrategy`: uniform over the enabled list.
* `BestFirstStrategy(GoalDistance)`: scores every enabled successor (or a
  random sample of them) by the goal distance after firing, computed with
  `Marking::peek` (apply, evaluate, revert), keeps the minimum with random
  tie-breaking, and takes an epsilon share of uniform random moves. Restarts
  after `stall` steps without improving the best distance of the run.
* `RelaxedPlanStrategy`: FF-style. Computes a relaxed plan from the current
  marking (below) and fires one of its helpful transitions; falls back to a
  random enabled transition when none is enabled. Same epsilon and stall
  mechanism, on the relaxed-plan value.

## Goal distances

* `MarkingDistance` (TAPAAL style): an atom's distance is how far its linear
  value is from satisfying the comparison; `And` sums, `Or` takes the minimum.
  Cheap, but blind to the structure: a token just past its target place looks
  close while it needs a full cycle.
* `StructuralDistance`: for atoms `place >= k` with no token on the place, the
  estimate becomes the number of transitions a token present in the marking
  needs to travel to the place, from a backward BFS over producers computed
  once per goal place. Optimistic: only one pre-place of each transition is
  accounted, so co-requirements (a control token that must be in a given
  state at the same time) are invisible.
* `RelaxedPlan` (planning h_add, as in directed unfolding): a Dijkstra pass
  from the marked places ignoring token consumption; a transition costs one
  plus the sum of its pre-place costs, a place takes the cheapest producer.
  The goal is evaluated on these costs (`And` sums, `Or` takes the cheapest
  branch, atoms other than `place >= k` fall back to the marking distance),
  and a relaxed plan is extracted backward through the recorded achievers.
  Plan transitions whose pre-places are all marked now are the *helpful*
  transitions. The pass stops as soon as every goal place has a cost; when
  the goal is far it visits the whole net, which on a net with 60k
  transitions costs a few milliseconds per step. Bookkeeping uses epoch
  stamps so nothing is cleared between steps.

* `BoundDistance`, for a bound target without a limit: a large offset minus
  the value of the form. Never zero (there is nothing to reach), monotone in
  the value, so best-first climbs the form and its stall detection and the
  shared pool see improvements. With a limit the bound is the atom
  `form >= limit` and every heuristic strategy applies unchanged; without a
  limit they all become this best-first.

* `DeadlockDistance`, for a deadlock target: the number of enabled
  transitions of the marking. Zero is exactly the target, so unlike
  `BoundDistance` this is a real distance and every best-first mechanism
  applies to it unchanged. It is the "fewest successors" rule of the 2020
  prototype, restated as a distance so it lives in the portfolio rather than
  beside it.

  Scoring a candidate `t` means counting the enabled transitions of
  `m + effect[t]` without building that marking. Only the consumers of the
  places `t` touches can change status, so the count moves by
  `-(disabled by t) + (enabled by t)` over that set, which is small on a
  sparse net: the same reasoning `EnabledSet::update` already applies after a
  real step, run here in peek mode over a sampled candidate set.

  A cheaper proxy, if the exact delta proves too slow: score `t` by how many
  transitions it disables, read off `consumers[p]` for the places it
  decreases, ignoring what it enables. That is the structural preference of
  `WALK_PLAN.md` 4.6 and needs no peek at all, but it is blind to a
  transition that disables three and enables four.

* `ParikhStrategy`: hint-driven. The net groups transitions by identical
  effect (`WalkNet::effectClassOf`, what the state equation sees); the hint's
  firing counts land on classes. A step chooses uniformly among the enabled
  transitions whose class still has count and decrements it; when none is
  enabled the run restarts. With probability one per mille per restart the
  restriction is skipped for the step, so a vector that is slightly off still
  guides. A focus with a hint runs the `--hintStrategies` pool.

The relaxed plan is what solves coordination-shaped goals (several
components that must each reach a specific state, gated by a shared control
token): the marking and structural distances see no gradient there.

## Walker: one run

```
restore initial marking and enabled snapshot; clear trace; check every open target
loop:
  if focus solved or no target open: stop
  if enabled empty: claim the deadlock targets; count a reset, restore, check all if pooled, continue
  if run length exhausted: count a reset, restore, check all if pooled, continue
  t = strategy.choose()  (RESTART: same as above)
  marking.apply(effect(t), enabled.onPlaceChanged, collect touched places); trace.push(t)
  check the open targets mentioning a touched place; claim those that hold
  every 64 steps: check the wall clock, the step budget, the stop flag
```

A claim calls back the owner (the portfolio) with the marking and, when
recording, the trace since the last reset; the callback prints the `FORMULA`
line at once. A strategy spec with the `+sat` suffix (`random+sat`, `bestfirst+sat:10:2000`)
puts the walker in saturation mode: the chosen transition is fired as many
times as the marking allows in one update (the largest k such that every
place it depletes keeps its input weight for k firings; 1 when it depletes
nothing). Stacks of tokens move in one step and places empty where a random
walk would rarely empty them; the trace records k firings and the strategy is
told through `onFired(t, k)` (the Parikh strategy consumes k). A poor choice
when only part of a stack should move, hence a modifier for some threads of a
portfolio, not the default.

Trace recording is optional (`WalkBudget::recordTrace`, CLI
`--trace`). When off, the walker never touches the trace vector and a witness
is reported as the marking reached. When on, the trace holds the transitions
fired since the last reset, so a witness is the sequence from the initial
marking; it is replayed on a fresh marking and checked before it is printed.
A witness found from a pooled restart has no trace.

## Deadlock

The target's `reached` test is already "the enabled list is empty", so a
deadlock target needs no new machinery in `TargetSet` or in the walker; it
needs a gradient. Without one the search is a uniform random walk, and on a
net whose dead markings sit off the bulk of the state space that walk does
not merely fail to conclude, it never comes near: 489 million steps and 489
resets on `Angiogenesis-PT-20` report zero dead ends.

With `DeadlockDistance` the deadlock target is an ordinary best-first target:
descend the enabled count, epsilon share of uniform moves, restart after
`stall` steps without improving the best count of the run. The epsilon and
the restart carry the weight here, because the rule is greedy and a net can
hold a region that is cheap in enabled transitions and has no dead marking in
it at all; `random` stays in the pool for the same reason.

Saturation belongs in any pool, not only this one, and it is free where it
cannot help: `Walker::maxFirings` returns 1 for a transition that depletes
nothing, and on a one-safe net the tokens beyond the first firing are zero, so
a saturated strategy degenerates to its plain form. The default property pool
saturates its three heuristics for that reason and keeps one plain random walk.

Saturation is the other axis of the deadlock pool, and on the evidence the
more decisive one. A net
whose dead markings sit behind a long repetition of one transition is drained
by firing that transition while it stays enabled, not by choosing it again and
again: on `Angiogenesis-PT-20` no unsaturated walk finds the deadlock in ten
seconds, and every saturated one finds it in milliseconds. The two axes are
independent, so the default pool of `--findDeadlock` spans both:
`random`, `random+sat`, `deadlock+sat:10:2000`, `deadlock:30:500`. On the three
models measured so far each of the first three wins one: `random+sat` takes
Angiogenesis in 5 ms, plain `random` takes Philosophers-PT-000100 and -001000.

A Parikh hint changes the picture entirely on a net whose deadlock is
structured rather than stumbled into, and `--findDeadlock` reads `--hints` for
that reason: the vector of firing counts that the caller's state equation
produces is the answer, not a nudge towards it. On `Philosophers-PT-005000` no
walk of any kind finds the deadlock in sixty seconds, and the hint "fire each
philosopher's first fork once" finds it in 39 ms, in exactly 5000 steps, the
minimal witness. There is one target in this mode and the caller names it after
its own property, so the hint is taken by count when the name does not match.
With a hint the pool becomes `parikh,deadlock+sat:10:2000,random+sat,random`.

The scale at which this stops working is worth knowing: the same hint on
`Philosophers-PT-010000` fires most of its counts and then restarts, 568 times
in a minute, without closing. No tool of the 2026 contest answered that
instance either.

`--findDeadlock` builds its strategies through `--strategies`, else
`--strategy`, else that pool, so a request for several threads spreads them
over it instead of running one strategy several times.

## Portfolio (threads)

`runPortfolio` starts N threads on one `TargetSet` and one focus; thread i
runs its own `Walker` with its own `Marking`, `EnabledSet`, RNG (seed + 7919 i,
the seed itself derived by the driver from `--seed`, default 1, per walk call:
two runs repeat up to thread timing, the shared claims, pool and counters)
and a fresh strategy instance built from `specs[i mod |specs|]` on the focus
goal (random when there is no focus), so all hot state is thread-local. The
`WalkNet`, the `TargetSet` and the goal expressions are shared read-only,
except the atomic solved flags. Claims are published under a mutex, in claim
order, through a callback; when the focus is claimed an atomic stop flag is
raised, which every walker polls every 64 steps (a step is milliseconds on a
hub-dense net, and the deadline overshot by seconds at 1024). After the threads join,
recorded traces are replayed by one verifier walker before being reported.
Strategy specs are `name[:epsilon[:stall]]`, e.g. `relaxed:0:300`.

The driver (`cli/WalkDriver.h`) applies a policy on top: an optional sweep
round (all threads random, no focus, all targets open) followed by rounds with
one focus per open target and a per-property budget growing tenfold per round
under `--totalTime`. Which mix works is an empirical question; the options
`--sweepTime`, `--roundTime` and `--strategies` are the knobs.

## Shared pool of restart states (optional)

`SharedPool` is a mutex-protected, bounded set of markings with their
heuristic value. At the end of a run (dead end, run length, stall) a walker
whose strategy exposes a best-of-run state (`Strategy::bestOfRun`) offers it;
the pool keeps it if there is room or it beats the worst entry, and refuses
duplicates. A restarting walker draws an entry with probability `shareProb`
instead of the initial marking, recomputing its enabled set from scratch. A
run that started from an entry reports whether it improved on the entry's
value; an entry that fails three times is evicted, and one a strategy finds
hopeless on arrival (`WalkContext::badStart`, the quest sweep with no target
questable from it) is evicted at once and the walk restarts again: the
defence against sharing states doomed by an earlier irreversible choice. A witness
found from a pooled restart has no trace (the pool does not keep traces), so
it is reported as a marking.

Random strategies neither publish nor benefit beyond drawing.

## Memory, rarity and restarts

`Knowledge` is what the threads know together: one relaxed atomic counter per
transition, the firings by anyone. `NoveltyTracker`, one per walker, holds the
walker's own counts, the transitions it fired at least once, the steps since it
last fired a new one, and the marking reached by the last *rare event*, a
transition nobody had fired before. Its counts go to the shared knowledge at
every reset and every 4096 steps, over the transitions touched since the last
merge only, so a small net pays nothing and a large one pays for what it did.

`RareStrategy` is the choice of a sweep: draw a few enabled transitions
uniformly (`--sample`, 8 by default) and fire the one fired the fewest times,
shared and own counts added, the longest enabled among equals (`EnabledSet`
stamps the clock at which a transition joined the list). O(sample) per step,
no scan of the enabled list. Uniform choice makes a transition that needs a
rare conjunction exponentially rare; this fires it the moment it is enabled.

`RestartPolicy` decides when a run ends, asked at every step with the run's
steps, wall clock and steps since novelty: a step budget, a wall clock budget,
a novelty stall, or any of them. A sweep restarts on all three (`--runLength`,
`--runTime`, `--noveltyStall`); a focused round keeps the step budget alone, a
deep hunt must not be cut by the clock; a strategy may still answer `RESTART`
on its own. A policy also says whether its restart *keeps the marking*: the
wall clock does, the strategy forgets (counters, tabu, current target) and
the walk goes on from where it stands, a reachable state and, for a quest
sweep, all of its progress; the step budget and the novelty stall restore the
initial marking or a pooled state, since leaving the state is their point. The marking of the last rare event goes
to the shared pool at the reset, so other threads restart from new ground.

Measured on ResIsolation-PT-N10P4, 1 000 `fireable` targets, 4 threads, 20 s:
restarts went from 0 to 19 per thread and the distinct transitions fired from
about 90 to 153 per thread, with the rare and the uniform choice alike. 153 is
where every 147 000-transition instance of the family stops, on the cluster as
here: a gate of the model (WALK_PLAN.md section 10), not of the choice.

## Components and quests

`Components` reads the net as sequential processes from its P-semiflows
(computed in-process by the invariant engine when a strategy asks,
`--strategy=sync`; flows with negative coefficients left places uncovered on
ErlangenMainframe, semiflows cover every place): a semiflow is a component,
its places the local states, its value at the initial marking the tokens it
carries. Not every semiflow is a process, and there can be exponentially
many: the decomposition takes a *covering*, the smallest semiflows first
while each covers a place not yet covered (DLCflexbar-PT-8b has 3 040
semiflows sharing hub places; every transition touching a hub belonged to
thousands of processes, and the tables took 32 GB), and it is refused when
the transitions' touches on the covering exceed `MAX_WORK`. The invariant
engine runs under a box of a quarter of the walk's total time in both
phases (a phase 1 stopped there yields no basis: FamilyReunion's 113 229
flows never end their phase 2), and the sweep takes what the box left of the
total. A transition
that takes from a component and gives back to it is a local move of that
process; one touching several components synchronises them (its *sync
degree*). Each component has a local graph, and `questDistances(c, p)` gives
every place of `c` its shortest path to `p` by the process's own moves, or,
where only synchronisations lead there, the distance over every move plus a
large offset: a process that can walk alone is driven first, one behind a
barrier still knows which way to lean. A degree profile signature groups
components that look alike.

`ComponentStrategy` (`sync`) walks toward a goal that is a conjunction of
single-place atoms, one quest per atom: `place >= k`, `place > k` and
`place == k` are all quested as the lower bound they imply (the claim itself
is checked on the exact predicate; ITS-Tools poses fireability on a safe net
as `p == 1`), and the process holding the place is driven to it. The distance of a marking is the sum over the unsatisfied quests
of the process's distance to its place; the choice is greedy on the distance
after firing among a few sampled enabled transitions (`--sample`, 16), with
an epsilon share of random moves that never undo a satisfied quest. A
satisfied quest *freezes* its process: a transition taking tokens from the
place is refused, so the processes arrive one by one and wait, and a product
of arrival chances becomes a sum.

When the goal places lie behind barriers, the strategy works in *stages*,
recursively: a stack whose root is the goal's quests. A stage whose every
open quest is behind a barrier has no move of its own left and pushes a
sub-stage at once; a stage with some quests still able to walk waits for a
stall, fifty steps in which no quest crossed its barrier (the far side of the
distance fluctuates with every local move and cannot be read for progress).
The sub-stage is the barrier every process can reach alone (a pre-place
already holding its tokens costs nothing whatever its component can do) whose
outcome brings the stalled stage's quests closest, the cheapest to align among
equals; its pre-places are the new quests; it is fired the moment it is
enabled, whatever it consumes, and its stage popped. A barrier whose crossing
did not lower its parent's distance is tabu at that level for the run; a stage
that cannot be staged further pops with its barrier tabu at the parent; the
root failing, or a quest whose tokens cannot reach its place, ends the run.
Depth is bounded. The relaxed plan was tried for the stage choice and
rejected: the delete relaxation lets a process be in two places at once, so it
never sees that a barrier moves processes away and it cycles between two
barriers.

The stage choice scans every barrier with per transition tables of the
(component, local place) pairs it consumes and produces and per component
distance matrices computed once, so a stage costs milliseconds on 147 000
barriers. Measured on 4 threads: ResIsolation-PT-N10P4 (147 855 transitions,
14 processes behind a chain of 13 forks, 147 456 fourteen-way barriers),
three barriers that random, best-first and relaxed-plan walks never fired,
nor 1800 s of every engine on the cluster, fire in about 95 steps, 2.7 to
7.5 s including the parse; ErlangenMainframeV2-PT-bP10C08 (79 094
transitions, 28 processes, 78 737 barriers of degree 12 or 13), three
transitions the campaign left open fire in 4.7 to 11.3 s.

`QuestSweep` (`--sweepChoice=sync`) is the sweep as a sequence of quests:
with no focus given, it ranks the open targets by the summed component
distance of their atoms from where the processes stand, takes the target of
its thread's rank so that four threads chase four targets, and runs a
`ComponentStrategy` toward it; when the target is claimed, the quest gives up
or the run restarts, it picks again. Every target satisfied on the way is
claimed by the walker, so a quest that reaches a shared pre-set claims every
transition on it. Above a few thousand open targets a random sample of them
is ranked, the nearest of a few thousand being near enough; the quest sweep
restarts ten times less often than the others, its progress (the processes
forked and placed) living in the marking. A marking from which no open target
can be quested, every one hopeless from where the processes stand or given
up on this run, restarts the walk (a pooled start is condemned, above); when
the restart lands no better the run goes on rarity without ranking again,
since a ranking at every step over a filler walk is the cost of a quest for
nothing. Within a quest the random share of
moves is sideways only, a move that neither undoes a process in place nor
worsens the summed distance by more than one, and a process falling behind a
barrier during a stage revises the stage at once. On 4 threads, 30 s, full
QuasiLivenessAll sets: ResIsolation-PT-N10P4 51 029 of 147 855 (the
campaign's 1800 s of every engine: 153), ErlangenMainframeV2-PT-bP10C08
1 570 of 79 094 (321), StigmergyElection-PT-10b 3 704 of 120 600 (797).

A quest `place >= k` gathers: the tokens still missing come from every
component holding the place, the distance is the sum of the k nearest tokens'
local paths (a token per unit of marking), and a token already on the place is
frozen there. A bound target on one place is a staircase in `QuestSweep`: the
quest is `place >= known + 1`, known being the shared record or the marking,
re-quested a step higher while it climbs; a quest that cannot start marks its
target unquestable for the run and the sweep moves on rarity, so the quest
sweep is never slower than the rarity sweep where quests do not apply.
