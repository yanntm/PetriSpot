# `walk/` — explicit exploration

Sparse, per-thread explicit walk toward a goal state. `algorithm.md` describes
the data structures and the step.

* `WalkNet.h` — compiled read-only view: effects, consumers per place sorted
  by weight, effect classes (transitions with identical effect).
* `Marking.h` — sparse marking with in-place `apply(effect, onChange)`.
* `EnabledSet.h` — enabled transitions by delta (unsatisfied-arc counters).
* `Target.h` — goal predicate (normal-form expression) or deadlock.
* `TargetSet.h` — the open targets of a run: atomic solved flags, the
  place-to-targets index, the deadlock targets.
* `Strategy.h` — transition choice interface (`RESTART` sentinel), `RandomStrategy`.
* `GoalDistance.h` — distance interface, `MarkingDistance` (expression distance).
* `StructuralDistance.h` — hop-based refinement for `place >= k` atoms.
* `BestFirstStrategy.h` — greedy on a `GoalDistance`, epsilon, stall restart.
* `RelaxedPlan.h`, `RelaxedPlanStrategy.h` — delete-relaxation cost (h_add),
  relaxed plan extraction, helpful-transition choice.
* `ParikhStrategy.h` — choice restricted to the effect classes with firings
  left in a Parikh hint, relaxed as restarts accumulate.
* `NetStats.h` — structural histograms (`--netStats`).
* `Walker.h` — restart loop over a `TargetSet` with an optional focus, budget,
  stop flag, incremental target checks, optional trace, pooled restarts.
* `Portfolio.h` — strategy specs and factory, N threads on one target set,
  claims published as they happen, trace verification.
* `SharedPool.h` — bounded shared pool of promising markings for restarts.

## Extending

* **A new strategy**: implement `Strategy<T>::choose` (return a transition of
  `ctx.enabled`, or `RESTART`), optionally `onReset` and `bestOfRun` (needed
  for the shared pool). Register a name in `makeStrategy` (`Portfolio.h`) so
  `--strategy` / `--strategies` accept it. Keep all mutable state inside the
  strategy: one instance is created per thread.
* **A new goal distance**: implement `GoalDistance<T>::of(marking)`; plug it
  into `BestFirstStrategy` through `makeStrategy`. Distances must be zero
  exactly when the goal holds.
* **Hot loop**: `Walker::run` -> `Strategy::choose` -> `Marking::apply` ->
  `EnabledSet::onPlaceChanged` -> `TargetSet::targetsOf` on touched places. Anything added there must stay proportional to
  the arcs touched by the fired transition; per-transition or per-place scans
  belong in the strategy's restart path or in one-time tables of `WalkNet`.
* **Instrumentation**: `WalkStats` (steps, resets, arc visits, flips) is printed
  per thread by the driver; `--debugSteps=n` traces the first n relaxed-plan
  decisions on stderr; `--netStats` prints the net's arity and fan-out profile.
* **Testing**: development models live outside the repo in `bench/models/`
  (extract MCC archives there); hand-written property files are in
  `Petri/test/props/`; the MCC harness is `~/git/MCC-drivers` with
  `BK_TOOL=petrispot` or `petrispotxred`. Keep any run under about 15 s per
  invocation and redirect long logs to `Petri/test/logs/`.
