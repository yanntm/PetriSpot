# `walk/` — explicit exploration

Sparse, per-thread explicit walk toward a goal state. `algorithm.md` describes
the data structures and the step.

* `WalkNet.h` — compiled read-only view: effects, consumers per place sorted
  by weight.
* `Marking.h` — sparse marking with in-place `apply(effect, onChange)`.
* `EnabledSet.h` — enabled transitions by delta (unsatisfied-arc counters).
* `Target.h` — goal predicate (normal-form expression) or deadlock.
* `Strategy.h` — transition choice interface, `RandomStrategy`.
* `Walker.h` — restart loop, budget, witness trace and its verification.
