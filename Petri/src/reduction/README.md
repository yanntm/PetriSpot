# `reduction/` — native structural reductions

Native header-only structural reduction, with an initial implemented subset.
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
| `Coordinator.h` | Cleanup/fallback fixed point, 11 rule statistics, optional pass-level trace policy |
| `rules/DuplicateTransition.h` | One duplicate-transition rule |
| `rules/ConstantPlace.h` | One constant-place rule |
| `rules/EmptySiphon.h` | One empty-siphon rule |
| `rules/ImplicitForkJoin.h` | One structural implicit-place rule |
| `rules/ReachabilityRelevance.h` | One property relevance rule |
| `rules/TrivialPost.h`, `rules/PreAgglo.h`, `rules/PostAgglo.h` | Separate agglomeration rules |
| `rules/DuplicatePlace.h`, `rules/NoEffect.h`, `rules/SinkPlace.h` | Separate cleanup rules |
| `Properties.h` | Adapter for `expr::Property`, support and checked rewriting |

Each rule gets its own file and named type; the list illustrates the layout,
not a fixed inventory. Shared helpers contain mechanics, not multiple rules.
All nine goal names are accepted; the reference's full rule coverage and exact
schedule are not implemented yet. Current temporal goals use strong local
rules only; LIVENESS retains dead-transition obligations. STATESPACE retains
all place coordinates and duplicate transitions, preserving raw counts without
metadata reconstruction. Local clear/replace operations maintain sparse
transposes without index shifting; publication compacts to a normal net.
Create files as their rules arrive, keeping each responsibility roughly below
500 lines. The kernel depends on `core/`, not parsers, CLI, walkers, SMT, or
libHSC's calculus. The property and PNET adapters sit outside that kernel.
Graphical rendering belongs in an IO adapter, outside the reduction kernel.

Current related services: [sparse substrate](../core/README.md),
[property AST](../expr/README.md), [PNET records](../io/PNET.md).

`--reduce` opts query analysis into the module. It runs after property parsing
and before walk/CTL/LP net compilation, using the union of query supports.
`--reductionMs` limits reduction work (default 15000); `--reductionNoAgglo`
disables agglomeration. `--trace`, input hints and LP hint export keep the
original net until lifting is implemented. Export and invariant-only requests
still describe the original net. The original input remains available.

Deferred: full reference schedule/rules, counting-record and image adapters,
application-level bounded visual capture/PDF, original witness lifting, and
libHSC vendoring. `NoTrace` compiles observation hooks out; enabled policies
currently observe whole rule passes, not individual applications. See the
implementation report [PS_REDUCTIONS.md](../../../PS_REDUCTIONS.md) and the
isolated [validation folder](../../test/reduction/README.md).
