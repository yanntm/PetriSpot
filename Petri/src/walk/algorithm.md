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

## Target

The goal predicate of the run: either a normal-form `Expression` evaluated on
the marking (`reached` is `goal.eval(marking)`), or deadlock (`reached` is
"the enabled list is empty"). Evaluation is currently from scratch at each
step; an incremental evaluation indexed by touched places is the next step
once profiling shows it matters.

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

The relaxed plan is what solves coordination-shaped goals (several
components that must each reach a specific state, gated by a shared control
token): the marking and structural distances see no gradient there.

## Walker: one run

```
restore initial marking and enabled snapshot; clear trace
loop:
  if target reached: return witness (trace, marking)
  if enabled empty or run length exhausted: count a reset, restore, continue
  t = strategy.choose()
  marking.apply(effect(t), enabled.onPlaceChanged); trace.push(t)
  every 1024 steps: check the wall clock and the step budget
```

The trace holds the transitions fired since the last reset, so a witness is
the sequence from the initial marking. Before it is reported the walker
replays it on a fresh marking and checks the goal again.
