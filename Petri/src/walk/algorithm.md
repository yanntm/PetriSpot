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
ones; `RandomStrategy` picks uniformly. The context exposes the marking, the
enabled set, the net and the RNG. Heuristic strategies plug in here.

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
