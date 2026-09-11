# Structural reduction design

Status: target design, with an initial subset implemented; the folder README
and `PS_REDUCTIONS.md` distinguish current behavior from planned capabilities.
The reference inspection
is in [itstools-review.md](itstools-review.md). This document defines a native
module intended to express ITS-Tools' reduction capabilities through named
configurations, while separating rules, sparse edits and coordination.
ITS-Tools is the profiled, extensively tested reference engine. Preserve its
special cases, fast paths and operational limits until measurements justify
changing them; unusual models and downstream metrics matter here.

## 1. Boundary and purpose

Input is an ordinary weighted P/T net: nonnegative initial marking, positive
stored arc weights, no inhibitor/reset arcs, and separate pre/post matrices.
An absent arc is zero. Incidence alone is insufficient: it loses read guards.
The result is another `SparsePetriNet<T>` usable by existing analysis engines.
No serialization round trip or second public net hierarchy is needed.

Java's dynamic-product `image`/`keepImage` mechanism is not in the first
implementation increment. Keep an optional extension point for that LTL
configuration rather than permanently excluding the capability. Ordinary
index remapping and name-based traceability are supported independently.
Executable witness lifting is a possible later capability.

The desired default integration is one reduction before building walk tables,
LP constraints, invariants, or libHSC events and variable order. Parse properties
against the original net first. Analysis caches belong to a particular net
revision and must be rebuilt after reduction. libHSC's NUPN tree must be
remapped for pure deletion or rebuilt after fusion; a fused place does not
automatically retain its original unit's safety guarantee.

Keep preprocessing single-threaded initially, with owned mutable state and no
global cache. Independent calls may run concurrently. Once published, a result
is immutable to its analysis consumers.

## 2. Proposed API and ownership

Illustrative signatures, not headers committed in advance of review:

```cpp
namespace petri::reduction {
template<class T> struct Request;
template<class T> struct Result;

template<class T>
Result<T> reduce(SparsePetriNet<T> net, const Request<T>& request);

template<class T>
PropertyResult<T> reduceForProperties(
    SparsePetriNet<T> net, std::span<const expr::Property> properties,
    const Options& options);
}
```

Passing an lvalue makes the intentional copy; passing `std::move(net)` transfers
ownership. The workspace takes that ownership without another full net copy.
The property overload derives the contract and returns rewritten properties
together with the net, so callers cannot accidentally mix coordinate systems.
It delegates to the same kernel as net-only analysis. No logging to stdout;
the caller formats optional progress and diagnostics.

`Request<T>` contains:

* A named reduction configuration, its preservation contract, and sparse set
  of observed original place IDs.
  No implicit reachability mode for an empty request: net-only callers must
  explicitly request their intended semantics.
* Required artifacts: selected counting measures (and, potentially later,
  original-trace replay), or
  neither. Optional artifacts may be invalidated with a recorded reason;
  required artifacts veto a rule that cannot maintain them.
* Known facts about this input revision, such as proven safety/bounds, and
  optional provenance records. Unknown is distinct from false. Caller-supplied
  facts are trusted premises, not inferred from a filename or NUPN label.
* Options/limits: enabled rules, schedule, agglomeration limits, search cutoffs,
  and optional deadline/cancellation. These control effort, never weaken
  semantic guards. No price or cost-model service is required.
* Optional reduction tracing, disabled by default, with a sink and explicit
  capture limits. This observes rule applications for study and visualization;
  it is separate from name-based traceability and executable witness lifting.
  See [tracing.md](tracing.md) for the capture and PDF design.

`Result<T>` contains:

* Owned reduced net, the effective preservation contract, and net revision.
* Original-to-result place mapping with explicit `kept(index)`, `constant(T)`,
  or `unavailable` entries; result-to-original identity/aggregate provenance.
  Deleted unobserved places need not be reconstructible. Names retain their
  traceability role; numeric maps serve sparse indexing and property rewriting.
* Updated counting records/facts, including reasons
  for invalidation. Renumbering never silently carries a stale record.
* Typed deductions (for example no reachable deadlock), with a rule and its
  premises. A dead transition proves non-liveness, not a reachable deadlock.
* Stop reason: fixed point of enabled/eligible searches, budget, cancellation,
  or deduction. Counts of attempts, matches, edits, removed/added places,
  transitions and arcs, plus elapsed time and skipped-work reasons.

Invalid input is an error before edits. Numeric overflow must not commit a
partial rewrite: reject the candidate with a diagnostic or terminate with the
last valid result and an explicit arithmetic-limit reason. The final choice
between these two behaviors can follow the existing arithmetic conventions.
Neither is a property verdict. A timeout returns a sound partial reduction.

## 3. Reduction goals and configurations

Keep all ITS-Tools reduction types as public named goals: `NONE`, `DEADLOCK`,
`REACHABILITY`, `SI_LTL`, `LTL`, `LIVENESS`, `STATESPACE`, `LI_LTL`, `SI_CTL`.
Each resolves to a configuration: preservation requirements, eligible rules,
coordinator schedule, default options/limits, and optional artifact handling.
The names remain recognizable to existing callers. Their implementation need
not be a switch repeated through every rule.

Configurations can share a common local-rule set and compose reusable phases.
They are not a simple hierarchy of increasingly strong reductions: liveness,
branching, trace length and counting impose different restrictions. A derived
configuration explicitly adds/removes rules or tightens their guards. Preserve
LI_LTL as a distinct goal with the reference's eligibility choices; do not
collapse it into SI_LTL based on a generic stuttering label.

The table below describes semantic building blocks, not replacements for the
named goals. Each rule documents a goal applicability table plus structural
guards. Goal support grows incrementally; a partially implemented configuration
reports its rule coverage, rather than silently aliasing a different goal.

| Profile | Required relation | Initial availability |
|---|---|---|
| Observed reachability | Equal sets of reachable valuations of observed places | First target |
| Deadlock existence | Original has a reachable deadlock iff result does | First target, separate rule guards |
| Strong behavior | Preserve branching, steps and deadlock under stated observation map | Conservative local subset for CTL with next |
| Stutter behavior | Appropriate divergence/deadlock-sensitive relation for the supported temporal fragment | Later, rule-by-rule justification |
| Counting | Exact requested original state/arc/token measures via maintained reconstruction | Local subset first |

Observed reachability preserves arbitrary state predicates over the retained
support, both EF and AG truth, and extrema of expressions on that support.
It does not license answering arbitrary CTL, transition liveness, or counting
the original states. A conjunction of obligations uses their intersection of
eligible rules; incompatible query families can instead use separate reduced
nets. Start with a shared union of supports, not one net per formula.

Protection means preserving the observed value, not preserving a numeric ID
or forbidding every adjacent arc edit. A proven constant can be substituted;
otherwise observed places remain explicit in the first version. A transition
is visible when its **effect** on the observation support is nonzero. A read
arc alone is not a visible update, although its guard remains semantically
essential. Do not confuse visibility with adjacency.

`Properties.h` collects support recursively through booleans and CTL state
atoms. Fireability is already expanded into original input-arc comparisons:
rewrite those comparisons, never reinterpret a removed transition ID.
Deadlock atoms and EX/AX require dedicated eligibility, not just place support.
Unsupported formulas get no inferred permissive profile. Initially use the
strong local subset for general CTL; do not enable SI_CTL based only on absence
of EX/AX until its divergence and finite-path semantics are established.

The existing AST stores coefficients/constants as `long long`, not `T`.
Kernel mappings and arithmetic stay in `T`; an adapter must check conversion
and arithmetic when substituting constants. If a rewrite does not fit, retain
that observed place or report the unsupported rewrite before publishing the
pair of net and properties. Do not silently narrow a 128-bit marking.

## 4. Mappings, evidence and metadata

An index permutation is enough for retained places. A constant substitution
is enough for a removed constant. A sum of original places describing a fused
place is a forward image; it cannot recover each original value. Agglomeration
can hide intermediate markings and change the initial marking. These are
different relationships and must not share one unqualified `image` field.

The following trace design is a deferred option, not part of the first API.
If added, trace support must constrain the schedule when requested.
For duplicates, map to a representative original transition. Preserve names
and composed names for human traceability, as ITS-Tools does. Record merges,
removals and renames against those names; numeric slot changes must not erase
that history. Expose the reference's excessive-name-length limit as an option;
if names are shortened, retain an optional old/new-name record when traceability
is requested. No provenance DAG is required for ordinary reduction.

If executable witness lifting is added, a simple agglomerate can additionally
store a recipe such as `sequence(h, repeat(f,k))`. Record any pre-fired
initial prefix. More general fusion can require marking-dependent routing;
skip it under replay requirements until that lifting algorithm exists.
Replay on the original net and check the original goal before publishing an
original witness. A macro expansion alone does not establish all rule guards.
Reversing Parikh hints across compositions is likewise not a permutation:
invalidate unsupported hints explicitly and let the existing engine run
without them. This does not invalidate the property itself.

Counting metadata follows [PNET](../io/PNET.md), including the distinction
between an absent block and an empty, present identity record. Keep an optional
typed record beside the net; convert to/from `PNETIO<T>::Blocks` in an adapter,
without requiring the reducer to include a file codec.

* Duplicate transitions merge `TMULT` weights; proven dead transitions only
  reindex them. Removing a possibly enabled transition needs a justified ghost
  record or invalidates arc counting.
* Constant deletion records the removed marking for token maxima/totals.
  An existing `PDROP` must be carried forward; absence is not a license to
  invent provenance for an already reduced imported net.
* Free SCC fusion can carry a `PCOEF` state weight, but that does not preserve
  arc counts. Weighted state counting and maximum tokens in an original place
  each need their own established relation, not an assumption from token sums.
* Unknown named metadata is retained only for an unchanged net; after edits,
  drop it unless its producer supplied a valid maintenance contract.
* Safety, bounds, invariant bases, decomposition and compiled transitions are
  distinct from these records. Revalidate, transform, or invalidate each.

Required count preservation skips incompatible rewrites. Optional records may
be lost. The reference configuration retains ITS-Tools' metadata-sensitive
choices, including the special STATESPACE route. Other configurations may make
these priorities explicit in the request. Later ghost support
must be implemented in producer and consumer before advertising arc recovery.

## 5. Mutation and sparse execution

Use the existing `SparseArray<T>`, `SparseBoolArray`, and `MatrixCol<T>`.
They already provide merge arithmetic, transposition, column edits, and batch
row deletion. The missing service is coordinated net mutation: current
`SparsePetriNet` exposes mutable matrices but const name/marking vectors and
has no atomic bulk replacement of all its fields. Its matrix constructor
assigns generated names. Initially a private workspace may copy labels and
use existing named builders at publication; measure that sparse O(arcs) cost
before adding a small move-based bulk constructor upstream. Never mutate only
matrices and leave dimensions, names, or maximum arc facts inconsistent.
`getMaxArcValue()` currently returns `int` despite its `T` field; do not use
that narrowing accessor as a correctness premise.

A rule recognizes against a consistent revision and prepares an edit:
removed objects with reasons, replacement columns, initial-marking changes,
mapping/fact updates, and touched neighborhoods. Arithmetic and artifact
eligibility are checked before commit. Conflicting candidates are rechecked
after an intervening mutation; do not assume all matches in a scan coexist.

Canonical matrices are transition columns. Lazily materialize their transposes
for place-oriented passes. The workspace offers local replacement, retirement,
append, and explicit compaction; rules choose the appropriate operation without
duplicating transpose/index maintenance. Preserve efficient reference paths:
avoid invalidating a whole transpose when only a few sparse entries change.

**Clear/retire, append, compact later** is the preferred transition-edit path.
Retiring a transition clears both columns and removes its entries from the
materialized pre/post transposes by traversing the old supports. It marks the
slot inactive without shifting any other transition index. Appending a new
transition adds its columns and sparse transpose entries at the end; transpose
row counts grow accordingly. Replacing a column visits its old/new support
union. These operations cost work proportional to touched arcs rather than
renumbering the entire matrix for each match.

An inactive cleared column is **not** an active source/no-effect transition.
Maintain an explicit active set, and make every rule, graph scan and deduction
respect it. A genuinely active transition with empty pre/post remains a real
transition. Counts and growth limits use active objects; a separate storage
limit monitors accumulated retired slots. Avoid reusing retired slots within
a pass, so queued IDs cannot silently acquire another meaning.

Trivial agglomeration remains its own fast rule even if a general rule matches
the same pattern. It can redirect a surviving transition's output and retire
the continuation, preserving matrix shape/indices and most adjacency data.
It changes arc contents, but avoids structural row/column deletion and a
general composition allocation. Do not route it through a heavyweight edit
representation or per-candidate full-net validation. Shared edit helpers can
commit a small prepared local change directly after its checks.

Compact at configured phase boundaries, when retired storage crosses a limit,
or before publishing the result. Compact with one old-to-new map, repairing
matrix keys, observations, names and metadata together. Place retirement needs
the same discipline and can initially be batched in place-oriented passes.
Rebuild transposes at compaction if cheaper than reindexing them; retain the
reference's immediate batch-deletion path where it is advantageous. The public
result has ordinary compact indices and no inactive slots.

Candidates carry revision information; recheck affected candidates after edits
and invalidate all slot-based candidates at compaction. O(P+T) bookkeeping at
preprocessing boundaries is acceptable; no dense work enters firing loops.

Hash immutable column pairs for duplicates and effects for dominance; compare
full sparse vectors to resolve collisions. Narrow candidates by degree and
adjacency before expensive comparisons. Dependency graph expansion costs
sum over transitions of |pre|*|post|, even with sparse storage. Use a budget or
an equivalent bipartite traversal; never drop dependency edges to fit a limit.

## 6. Rules and scheduling

**One rule, one file, one named type.** For example,
`rules/DuplicateTransition.h` defines `DuplicateTransition`, and
`rules/TrivialPost.h` defines `TrivialPost`. A family such as place cleanup is
a scheduler grouping, never a class containing all its constituent rules.
Shared graph and composition helpers do not own rule selection or policy.
Each rule's file states its contract next to its implementation; introduce a
`rules/README.md` inventory when implementation creates that subfolder.

The proposed conceptual protocol is:

```cpp
bool eligible(const Contract& contract) const;
// Structural search yields a typed Match, or no match, under a work budget.
SearchResult<Match> find(const NetView<T>& net, const RuleOptions& options) const;
Edit<T> prepare(const NetView<T>& net, const Match& match) const;
```

These signatures express responsibilities; template placement and the exact
search cursor/result types remain design choices. No inheritance hierarchy or
runtime plugin registry is needed initially: named rule types and explicit
calls in the coordinator suffice. `find` distinguishes exhaustion from a
limit stop. Rules may implement a batch pass using these responsibilities
internally, so trivial rewrites need not allocate a generic match/edit object
per application. A match carries its source revision;
`prepare` rechecks validity if that revision changed.

The guard has two parts: eligibility for what the caller wants to preserve,
then structural premises on this candidate. Artifact maintenance is checked
before committing the edit. Both parts are mandatory even if the scheduler
already filtered the rule family. A rule may expose distinct named variants
when their proof obligations differ, rather than nested boolean switches.

Use explicit options and limits rather than a mandatory cost estimate. A rule
can cheaply decide to skip a search on the current net or reject an expansion
candidate: transition count, degree-bucket size, search depth, producer/consumer
product, generated arc count, applications per pass. The coordinator can also
inspect recent progress to decide whether to enter a phase. That is already
dynamic selection without assigning a speculative price to every rule.
Keep measured time and actual changes in statistics for later tuning.

Preserve the reference defaults and exact comparisons in its configuration:
implicit depth 5; composition skipped above 20,000 transitions; future buckets
skipped at 10,000; complex post cross-product rejected at 32 when both sides
are non-singleton; complex post stops after more than 100 applications; global
loop stops after four consecutive transition-growing rounds. The image-mode
consumer limit belongs to its optional configuration. Give these limits names,
units and meanings; changing one must not change a structural guard.

Each rule states: mathematical preconditions; supported profiles; observation
and artifact obligations; sparse recognition algorithm; edit; and validation
examples. Split cheap rejection tests from the semantic predicate and split
that predicate from construction. Share agglomerate construction, not an
opaque nest of `doComplex`/`doSimple` switches across different proofs.

The **rule coordinator** assembles and runs the selected goal's rules. Its
small vocabulary is ordered sequence, repeat while changed, fallback when a
phase made no progress, and conditional phase. Ordinary typed C++ functions
are sufficient; no scheduling language is required. The reference configuration
reproduces the order and nested loops in the review, including initial SCC
passes, trivial-before-general post, late siphons/token movement, and the
separate STATESPACE path. It is the initial comparison baseline. A modified
schedule is a named configuration, so experimental ordering does not silently
replace that baseline.

Each rule returns pass progress, limit/skipping reasons and statistics. The
coordinator owns fixed-point decisions and growth tracking, not individual
rules. It can rerun cheap cleanup after broader rewrites according to the
chosen schedule. Reusable phases preserve the reference's fast fallbacks.
Rules also expose optional before/after capture points with their local focus
and highlighted objects. Dispatch once into an untraced or traced coordinator;
the untraced specialization compiles out capture, descriptions and snapshots.
Tracing never requires keeping retired columns alive after a rewrite.
Use actual committed change, including guard or marking changes, as progress;
net size alone is not sufficient. End at stability only if every enabled
eligible search completed without a match; report exhausted candidate budgets
as limited search, not global irreducibility.

The reference schedule aims at operational compatibility, not byte-identical
Java hash iteration. Document candidate traversal and tie breaks; preserve
reference ordering where observable, and record differences. Keep thresholds
in named option groups, including optional arc/storage limits absent from the
reference. Such extra limits are disabled in the compatibility configuration
unless requested. Check cancellation between batches and before large allocation.

For example, post-agglomeration through initially empty p with unit arcs,
one producer h and one consumer f, h's only output p and f's only input p,
and both invisible, replaces h,f by h;f. Outside p its pre is pre(h), post is
post(f). Preconditions exclude h=f. The intermediate marking is hidden, so
this is not a next-step preservation rule. Broader weighted composition uses
k=post(h,p)/pre(f,p) only after divisibility and continuation conditions are
proved; build pre(h)+k*pre(f) and post(h)+k*post(f), omitting p, using checked
arithmetic. This formula is not a generic composition of arbitrary transitions.

## 7. Implementation stages after design review

1. Introduce named goals/configurations, coordinator, request/result, sparse
   local edits, retirement and compaction, names and remapping; implement
   duplicate transitions, proven constants and dead transitions. Preserve
   observations and basic counting records, with explicit refusal of any
   unsupported required artifact. Integrate one preprocessing point in
   PetriSpot, then vendor through libHSC's existing upstream-copy workflow.
2. Add empty siphons, local implicit-place detection and reachability relevance.
   Establish deadlock guards independently, especially source transitions.
3. Add the trivial chain above and simple pre/post agglomeration; compare
   schedule costs before expanding rule coverage. Decide separately whether
   original witness lifting warrants implementation.
4. Complete the structural rule inventory and goal-specific eligibility,
   including free SCC/future fusion and weighted/partial agglomerations.
   Bring each named goal's reference configuration to parity; document the
   original conditions and their rationale alongside each rule. Keep the
   optional image configuration as a separately scheduled capability.
5. Add optional invariant/LP/SMT-derived facts through a separate adapter if
   needed to reproduce the outer ITS-Tools orchestration. Structural
   preprocessing remains usable with no solver dependency. Explicit transforms
   such as read abstraction and causal decomposition remain separately callable
   operations with their own contracts, rather than default equivalence rules.

Maintain a compatibility checklist as implementation starts: each reference
rule/variant, goal gate, option, fast path, caller-side preprocessing requirement,
and metadata/name behavior is implemented, deferred, or intentionally changed
with a reason. Incremental delivery is not permanent loss of reference scope.
Include ITS-Tools' high-debug views in that checklist. Add optional capture
hooks alongside rules, then local DOT export and a multipage PDF renderer as a
separate feature increment. Bounded capture and disabled-path overhead checks
are acceptance criteria for tracing, not optional later optimizations.

Validation should target semantics rather than reproduce implementation.
For small bounded examples exhaust original and reduced reachable markings:
compare projected valuation sets, deadlock existence, or weighted counts as
the contract requires. Test near misses: read arcs, nonzero initial p,
indivisible weights, two competing consumers, visible intermediate states,
source transitions, duplicate chains, and overflow in composition. Test mixed
property supports, deletion followed by fusion, metadata reindexing and replay.
For CTL use formulas distinguishing branching, next steps, and deadlock.

Use ITS-Tools as the performance and functional baseline, alongside small exact
semantic checks. Compare the reference configuration first; investigate changed
reductions rather than dismissing different sizes as expected. Distinguish
candidate-order differences from missed rules and avoidable performance loss.
Measure reduction time, peak storage, arcs and downstream analysis cost on
user-selected MCC models. Retain pathological cases motivating particular
limits/fast paths as regression inputs or external corpus references, including
cases where a smaller net makes a downstream metric worse. Rule-level timings,
transpose rebuilds, compaction work and allocated columns help explain losses.
Individual diagnostics stay within the repository's 15-second bound. No tests
or benchmarks are needed for this documentation-only stage.

## 8. Review decisions

Keep all ITS-Tools goals as configurations, sharing rule sets and phases where
appropriate. The coordinator owns ordering, fallback and fixed points; each
rule owns its guards and efficient rewrite in one file/type. Named options and
limits control dynamic admission; there is no mandatory cost model. Preserve
trivial fast paths, name-based traceability, and clear/append edits with local
transpose maintenance and deferred compaction. The reference configuration is
the baseline; changes need evidence, including on difficult models. Image and
executable witness lifting remain deferred optional capabilities. libHSC reuses
the native API through vendoring; PNET stays an adapter. Implementation proceeds
incrementally under the approved design; the report records remaining coverage.

## 9. Implemented reachability/deadlock schedule

Each named rule owns its own recognition and edits. The coordinator retains
Java's nested phases: free SCC and prefix first; place cleanup, transition
cleanup, prefix, implicit places, conditional free SCC, trivial post or simple
post to stability. Then try simple pre, future equivalence, non-complex post,
complex post, complex pre, free SCC and prefix, with no-progress gates. Always
try siphon-enabled cleanup. At stability try redundant composition, simple free,
complex free, partial free, partial post (reachability), then initial token
movement and cleanup. Four consecutive transition-growing outer rounds stop.

The source-transition deadlock deduction is checked before graph analysis as
well as in transition cleanup, making its caller prerequisite explicit. Pure
reachability and deadlock are the completed rule-body scope; SI-specific token
movement/image cases, counting modes and external SMT rules are separate work.
Exact edit counts and SCC representatives are not semantic results: stable slots
and checked arithmetic replace Java's index-shifting mutations. Rule admissions,
phase gates and limits retain their reference meaning.

Graph algorithms and their precise pruning behavior are documented in
[graph/README.md](graph/README.md); each rule is mapped in
[rules/README.md](rules/README.md). Bounds use the reachability configuration
with the objective's places protected. The specialized positive-effect dominance
method remains unscheduled, matching its commented-out Java caller.
