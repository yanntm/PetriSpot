# `walk/` — explicit exploration

Sparse, per-thread explicit walk toward a goal state. `algorithm.md` describes
the data structures and the step.

* `WalkNet.h` — compiled read-only view: effects, consumers per place sorted
  by weight.
* `Marking.h` — sparse marking with in-place `apply(effect, onChange)`.
* `EnabledSet.h` — enabled transitions by delta (unsatisfied-arc counters).
* `Target.h` — goal predicate (normal-form expression) or deadlock.
* `Strategy.h` — transition choice interface (`RESTART` sentinel), `RandomStrategy`.
* `GoalDistance.h` — distance interface, `MarkingDistance` (expression distance).
* `StructuralDistance.h` — hop-based refinement for `place >= k` atoms.
* `BestFirstStrategy.h` — greedy on a `GoalDistance`, epsilon, stall restart.
* `RelaxedPlan.h`, `RelaxedPlanStrategy.h` — delete-relaxation cost (h_add),
  relaxed plan extraction, helpful-transition choice.
* `NetStats.h` — structural histograms (`--netStats`).
* `Walker.h` — restart loop, budget, stop flag, optional trace, pooled restarts.
* `Portfolio.h` — strategy specs and factory, N racing threads, first winner.
* `SharedPool.h` — bounded shared pool of promising markings for restarts.
