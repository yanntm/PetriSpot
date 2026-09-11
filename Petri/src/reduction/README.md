# `reduction/` — native structural reductions

Design proposal for review; there is no implementation in this folder yet.
The module will transform a sparse P/T net before analysis, with an explicit
contract describing which observations and answers survive the transformation.
PetriSpot owns the implementation; libHSC can vendor it alongside `core/`.

Read [algorithm.md](algorithm.md) for the proposed interface, correctness
boundary, mutation strategy, and implementation stages. Read
[itstools-review.md](itstools-review.md) for the reference inventory, actual
rule schedule, optimizations, and points requiring further audit.
Read [tracing.md](tracing.md) for optional high-debug local views and the future
multipage PDF showing reduction steps, with bounded capture and no trace work
in the disabled hot path.

Proposed source responsibilities (names are provisional):

| File | Responsibility |
|---|---|
| `Reduce.h` | Typed request/result and entry point; no rule bodies |
| `Workspace.h` | Owned sparse net, adjacency views, validated edits and compaction |
| `Contract.h` | Preservation requirements, observations, rule eligibility |
| `Configuration.h` | All ITS-Tools named goals, shared phases and options/limits |
| `Mapping.h` | Names and traceability, index remapping, constant substitutions |
| `Metadata.h` | Maintenance or invalidation of counting records and structural facts |
| `Coordinator.h` | Configured rule ordering, fallback, fixed points and statistics |
| `Trace.h` | Optional bounded observation of rule applications; disabled specialization |
| `rules/DuplicateTransition.h` | One duplicate-transition rule |
| `rules/ConstantPlace.h` | One constant-place rule |
| `rules/EmptySiphon.h` | One empty-siphon rule |
| `rules/ImplicitForkJoin.h` | One structural implicit-place rule |
| `Graph.h` | Sparse dependency traversal, SCCs and stabilization |
| `rules/ReachabilityRelevance.h` | One property relevance rule |
| `rules/FreeScc.h`, `rules/FutureEquivalent.h` | Separate fusion rules |
| `Compose.h` | Checked construction of agglomerated transitions |
| `rules/TrivialPost.h`, `rules/PreAgglo.h`, `rules/PostAgglo.h` | Separate agglomeration rules |
| `Properties.h` | Adapter for `expr::Property`, support and checked rewriting |

Each rule gets its own file and named type; the list illustrates the layout,
not a fixed inventory. Shared helpers contain mechanics, not multiple rules.
The reference configuration preserves ITS-Tools' schedule, goal distinctions,
limits and fast paths. Local clear/replace/append operations maintain sparse
transposes without repeated index shifting; explicit compaction publishes a
normal net. Alternative schedules are separate configurations.
Create files as their rules arrive, keeping each responsibility roughly below
500 lines. The kernel depends on `core/`, not parsers, CLI, walkers, SMT, or
libHSC's calculus. The property and PNET adapters sit outside that kernel.
Graphical rendering belongs in an IO adapter, outside the reduction kernel.

Current related services: [sparse substrate](../core/README.md),
[property AST](../expr/README.md), [PNET records](../io/PNET.md).
