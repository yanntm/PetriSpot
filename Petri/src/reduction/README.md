# `reduction/` — native structural reductions

Native header-only structural reduction, with the reachability/deadlock structural rule inventory implemented.
The module transforms a sparse P/T net before analysis, with an explicit
contract describing which observations and answers survive the transformation.
PetriSpot owns the implementation; libHSC can vendor it alongside `core/`.

Read [algorithm.md](algorithm.md) for the proposed interface, correctness
boundary, mutation strategy, and implementation stages. Read
[itstools-review.md](itstools-review.md) for the reference inventory, actual
rule schedule, optimizations, and points requiring further audit.
Read [tracing.md](tracing.md) for optional high-debug local views and the future
multipage PDF showing reduction steps, with bounded capture and no trace work
in the disabled hot path.

Implemented source responsibilities:

| File | Responsibility |
|---|---|
| `Reduce.h` | `reduce(net, configuration, support[, trace])`, result and index maps |
| `Workspace.h` | Owned sparse net, adjacency views, validated edits and compaction |
| `Configuration.h` | All ITS-Tools named goals, options/limits and agglomeration eligibility |
| `Counting.h` | The counting record (transition multiplicities, place coefficients, dropped constants) and its maintenance under STATESPACE |
| `Coordinator.h` | Java phase order, 26 rule/phase statistics, optional pass-level trace policy |
| `rules/DuplicateTransition.h` | One duplicate-transition rule |
| `rules/ConstantPlace.h` | One constant-place rule |
| `rules/EmptySiphon.h` | One empty-siphon rule |
| `rules/ImplicitForkJoin.h` | One structural implicit-place rule |
| `rules/ReachabilityRelevance.h` | One property relevance rule |
| `rules/TrivialPost.h`, `rules/PreAgglo.h`, `rules/PostAgglo.h` | Separate agglomeration rules |
| `rules/DuplicatePlace.h`, `rules/NoEffect.h`, `rules/SinkPlace.h` | Separate cleanup rules |
| `Properties.h` | Adapter for `expr::Property`, support and checked rewriting; one reduction step for a property set |
| `Pipeline.h` | `prepare(...)`: constants into the formulas, initial state, decided properties dropped, reduce, to a fixpoint |
| `PropertyFacts.h` | Constant substitution and formula simplification |
| `Composition.h`, `TransitionAlgebra.h` | Checked sparse composition and effect operations |
| `graph/` | SCCs, dependency prefixes and stabilizing analysis |
| `rules/FreeSCC.h`, `rules/PrefixOfInterest.h`, `rules/LoopBack.h` | Separate graph rules |
| `rules/FreeAgglo.h`, `rules/PartialFreeAgglo.h`, `rules/PartialPostAgglo.h` | Full and partial free/post rules |
| `rules/FutureEquivalent.h`, `rules/RedundantComposition.h` | Future fusion and transition dominance/composition |
| `rules/ScalarTransition.h`, `rules/SinkTransition.h`, `rules/SourceTransition.h` | Transition cleanup and deadlock deduction |
| `rules/InitialTokenMove.h` | Final-stability pre-firing |
| `rules/DeadTransition.h` | Transitions the state equation proves never enabled (`lp/DeadTransitions.h`), budgeted by `deadMs` |
| `rules/BoundsDominance.h` | Unscheduled Java bounds-specific method |
| `cli/` | Standalone model/formula transformation, the counting record as PNET blocks, solved-property reporting |

Each rule gets its own file and named type; the list illustrates the layout,
not a fixed inventory. Shared helpers contain mechanics, not multiple rules.
All nine goal names are accepted. Reachability/deadlock follow the Java
structural phases; full SI-mode validation and caller-level SMT orchestration
remain outside the completed scope. LIVENESS retains dead-transition obligations. STATESPACE
runs the rules audited for the counting record (constant places, duplicate
transitions, no-effect transitions once arcs are untracked, free SCC, dead
transitions by the state equation) and the
workspace maintains the record through them: `TMULT` while the arcs are those of
the input, `PDROP` for removed constant places, `PCOEF` for fused free
components (algorithm.md section 4, `io/PNET.md`).
Constant representatives of fused components retain their coefficient as
isolated places; recording only their token sum in `PDROP` would lose states.
Local clear/replace operations maintain sparse transposes without index
shifting; publication compacts to a normal net.
Create files as their rules arrive, keeping each responsibility roughly below
500 lines. The kernel depends on `core/`, not parsers, CLI, walkers, SMT, or
libHSC's calculus. The property and PNET adapters sit outside that kernel.
Graphical rendering belongs in an IO adapter, outside the reduction kernel.

Current related services: [sparse substrate](../core/README.md),
[property AST](../expr/README.md), [PNET records](../io/PNET.md).

Every engine starts with `Pipeline.h`'s `prepare`, `--reduce` or not: the
constant places of the net (`PropertyFacts.h`: unchanged by every transition, or
in the greatest initially empty siphon) are substituted into the formulas, the
formulas simplified and confronted with the initial marking
(`expr/InitialState.h`: decided formulas answered, an until whose left side
fails initially replaced by its right side, `EF p` / `AG p` requalified as
reachability kinds), and the decided properties reported and dropped. With
`--reduce` the net is then reduced for the kinds and the union support that
remain, and the round repeats on the reduced net while something changes: a
smaller support lets more rules fire, a reduced net exposes more constants.
Without it only the formula side runs, the degraded mode. `--reduce` runs before
walk/CTL/LP net compilation.
`--reductionMs` limits reduction work (default 15000); `--reductionNoAgglo`
disables agglomeration; `--deadMs` (default 3000, 0 disables) budgets the
state-equation dead transition tests, which open every outer round and the
STATESPACE loop, and are reported as one diagnostics line (found, tested,
solves, pivots, passes, ms, budget hit). A pass cut by the budget resumes
after the last transition tested. `--trace`, input hints and LP hint export keep the
original net until lifting is implemented. Export and invariant-only requests
still describe the original net. The original input remains available.

Deferred: SI-mode completion, caller-level SMT orchestration, image adapters,
application-level bounded visual capture/PDF, original witness lifting, and
libHSC vendoring. `NoTrace` compiles observation hooks out; enabled policies
currently observe whole rule passes, not individual applications. See the
implementation report [PS_REDUCTIONS.md](../../../PS_REDUCTIONS.md) and the
isolated [validation folder](../../test/reduction/README.md).

Standalone transformation: `petri64 reduce -i MODEL --props FORMULAS --output DIR`
writes a matched PNET/formula pair plus names; `petri64 reduce -i MODEL --goal
STATESPACE --output DIR` writes the net reduced for counting with its record as
PNET blocks, for `hsc-pn --net DIR/model.pnet --states`; see [cli/README.md](cli/README.md).
Normal MCC analysis reports and drops constant properties before creating goals.
Structural deadlock deductions are explicit results, independent of net mutation.
