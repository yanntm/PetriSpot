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
* `targetsOf(p)`: the targets whose goal mentions place `p`, built once
  (sparse: only the places that appear in some atom have a non-empty list);
* `deadlocks()`: the deadlock targets, checked only when the enabled list is
  empty;
* per bound target, an atomic running maximum: a walker publishes a value
  only when it beats the last one it published (a thread-local array), so
  the atomic is touched at improvements only; the driver prints it at exit.

A walker re-evaluates a target only when the fired transition changed a place
the goal mentions: the places touched by `Marking::apply` are collected, their
target lists merged with an epoch stamp per target (no clearing), and each
candidate still open is evaluated in full (goal expressions are small trees).
The cost per step is proportional to the arcs touched and to the targets that
depend on them, never to the number of targets. Every open target is evaluated
in full once at the start of a run and after each reset (the marking is then
the initial one or a pooled state).

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
  if enabled empty: claim the deadlock targets; count a reset, restore, check all, continue
  if run length exhausted: count a reset, restore, check all, continue
  t = strategy.choose()  (RESTART: same as above)
  marking.apply(effect(t), enabled.onPlaceChanged, collect touched places); trace.push(t)
  check the open targets mentioning a touched place; claim those that hold
  every 1024 steps: check the wall clock, the step budget, the stop flag
```

A claim calls back the owner (the portfolio) with the marking and, when
recording, the trace since the last reset; the callback prints the `FORMULA`
line at once. Trace recording is optional (`WalkBudget::recordTrace`, CLI
`--trace`). When off, the walker never touches the trace vector and a witness
is reported as the marking reached. When on, the trace holds the transitions
fired since the last reset, so a witness is the sequence from the initial
marking; it is replayed on a fresh marking and checked before it is printed.
A witness found from a pooled restart has no trace.

## Portfolio (threads)

`runPortfolio` starts N threads on one `TargetSet` and one focus; thread i
runs its own `Walker` with its own `Marking`, `EnabledSet`, RNG (seed + 7919 i)
and a fresh strategy instance built from `specs[i mod |specs|]` on the focus
goal (random when there is no focus), so all hot state is thread-local. The
`WalkNet`, the `TargetSet` and the goal expressions are shared read-only,
except the atomic solved flags. Claims are published under a mutex, in claim
order, through a callback; when the focus is claimed an atomic stop flag is
raised, which every walker polls every 1024 steps. After the threads join,
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
value; an entry that fails three times is evicted, which is the defence
against sharing states doomed by an earlier irreversible choice. A witness
found from a pooled restart has no trace (the pool does not keep traces), so
it is reported as a marking.

Random strategies neither publish nor benefit beyond drawing.
