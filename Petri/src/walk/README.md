# `walk/` — explicit exploration

Sparse, per-thread explicit walk toward a goal state. `algorithm.md` describes
the data structures and the step.

* `WalkNet.h` — compiled read-only view: effects, consumers per place sorted
  by weight, effect classes (transitions with identical effect).
* `Marking.h` — sparse marking with in-place `apply(effect, onChange)`.
* `EnabledSet.h` — enabled transitions by delta (unsatisfied-arc counters).
* `Target.h` — goal predicate (normal-form expression) or deadlock.
* `TargetSet.h` — the open targets of a run: atomic solved flags, the places
  each target mentions with the direction that can make it hold, the
  deadlock targets.
* `Strategy.h` — transition choice interface (`RESTART` sentinel), `RandomStrategy`;
  the `WalkContext` carries the shared `Knowledge` and the walker's own counts.
* `RareStrategy.h` — the sweep's choice: the least fired of a few sampled
  enabled transitions, the longest enabled among equals.
* `Knowledge.h` — shared firing counts, one relaxed atomic per transition.
* `NoveltyTracker.h` — a walker's memory: own counts merged into `Knowledge`,
  coverage, steps since a new transition, the marking of the last rare event.
* `RestartPolicy.h` — when a run ends: step budget, wall clock, novelty stall,
  any of them; injected into the walker.
* `Components.h` — the net as processes, from its P-flows: local graphs,
  local distances to a place, synchronisation degree of transitions, a
  signature for look-alike components.
* `ComponentStrategy.h` — `sync`: one quest per goal atom, processes driven
  to their places and frozen there, stages through the barriers whose outcome
  brings the goal closest, a tabu on the barriers that did not help.
* `QuestSweep.h` — the sweep as quests: each thread picks the open target of
  its rank by component distance and runs `sync` toward it, again and again.
* `GoalDistance.h` — distance interface, `MarkingDistance` (expression distance).
* `StructuralDistance.h` — hop-based refinement for `place >= k` atoms.
* `DeadlockStrategy.h` — greedy descent of the enabled-transition count, the
  distance of a deadlock target; successor counts from the unsatisfied-arc
  counters, without building the successor marking.
* `BestFirstStrategy.h` — greedy on a `GoalDistance`, epsilon, stall restart.
* `RelaxedPlan.h`, `RelaxedPlanStrategy.h` — delete-relaxation cost (h_add),
  relaxed plan extraction, helpful-transition choice.
* `ParikhStrategy.h` — choice restricted to the effect classes with firings
  left in a Parikh hint, relaxed as restarts accumulate.
* `NetStats.h` — structural histograms (`--netStats`).
* `Walker.h` — restart loop over a `TargetSet` with an optional focus, budget,
  stop flag, incremental target checks over a thread-local up/down index of
  the targets it owns, optional trace, pooled restarts.
* `Task.h` — an exploration task for the scheduler: slice, slice report, share
  and virtual clock, totals.
* `WalkTask.h` — a task over one `Walker` and its strategy, resumable; the
  `ThreadReport` it closes into.
* `Scheduler.h` — time sharing: a heap of runnable tasks on their virtual
  clock, runner threads taking one slice each (steps, capped by the clock).
* `Portfolio.h` — strategy specs and factory, the tasks built from them and
  run by the scheduler on one target set, claims published as they happen,
  the per-kind summary, trace verification.
* `SharedPool.h` — bounded shared pool of promising markings for restarts.

## Extending

* **A new strategy**: implement `Strategy<T>::choose` (return a transition of
  `ctx.enabled`, or `RESTART`), optionally `onReset` and `bestOfRun` (needed
  for the shared pool). Register a name in `makeStrategy` (`Portfolio.h`) so
  `--strategy` / `--strategies` accept it. Keep all mutable state inside the
  strategy: one instance is created per task; what tasks share goes
  through `Knowledge`, read from the `WalkContext`.
* **A new restart rule**: implement `RestartPolicy::shouldRestart` over the
  `RunView` and add it to the `AnyOf` the driver builds (`WalkDriver.h`).
* **A new goal distance**: implement `GoalDistance<T>::of(marking)`; plug it
  into `BestFirstStrategy` through `makeStrategy`. Distances must be zero
  exactly when the goal holds.
* **Hot loop**: `Walker::run` -> `Strategy::choose` -> `Marking::apply` ->
  `EnabledSet::onPlaceChanged` -> the walker's own `up`/`down` lists of the
  touched places. Anything added there must stay proportional to the arcs
  touched by the fired transition; per-transition or per-place scans belong
  in the strategy's restart path or in one-time tables of `WalkNet`.
* **Instrumentation**: `WalkStats` (steps, resets, arc visits, flips) is printed
  per thread by the driver; `--debugSteps=n` traces the first n relaxed-plan
  decisions on stderr; `--netStats` prints the net's arity and fan-out profile.
* **Testing**: development models live outside the repo in `bench/models/`
  (extract MCC archives there); hand-written property files are in
  `Petri/test/props/`; the MCC harness is `~/git/MCC-drivers` with
  `BK_TOOL=petrispot` or `petrispotxred`. Keep any run under about 15 s per
  invocation and redirect long logs to `Petri/test/logs/`.
