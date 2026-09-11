# Structural reduction design

Status: proposal for review, before implementation. The reference inspection
is in [itstools-review.md](itstools-review.md). This document defines a native
module, not a promise to reproduce every Java transformation or its schedule.

## 1. Boundary and purpose

Input is an ordinary weighted P/T net: nonnegative initial marking, positive
stored arc weights, no inhibitor/reset arcs, and separate pre/post matrices.
An absent arc is zero. Incidence alone is insufficient: it loses read guards.
The result is another `SparsePetriNet<T>` usable by existing analysis engines.
No serialization round trip or second public net hierarchy is needed.

Java's dynamic-product `image`/`keepImage` mechanism is outside the proposed
scope. Its specialized LTL use does not justify carrying it into the native
interface. Ordinary index remapping is necessary independently. Original
witness lifting is a possible later capability, not an initial requirement.

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
ownership. Do not additionally copy the entire net into a Java-like reducer.
The property overload derives the contract and returns rewritten properties
together with the net, so callers cannot accidentally mix coordinate systems.
It delegates to the same kernel as net-only analysis. No logging to stdout;
the caller formats optional progress and diagnostics.

`Request<T>` contains:

* A preservation contract and sparse set of observed original place IDs.
  No implicit reachability mode for an empty request: net-only callers must
  explicitly request their intended semantics.
* Required artifacts: selected counting measures (and, potentially later,
  original-trace replay), or
  neither. Optional artifacts may be invalidated with a recorded reason;
  required artifacts veto a rule that cannot maintain them.
* Known facts about this input revision, such as proven safety/bounds, and
  optional provenance records. Unknown is distinct from false. Caller-supplied
  facts are trusted premises, not inferred from a filename or NUPN label.
* Options: enabled rule families, schedule, deadline/cancellation, work and
  expansion budgets. These control effort, never weaken semantic guards.

`Result<T>` contains:

* Owned reduced net, the effective preservation contract, and net revision.
* Original-to-result place mapping with explicit `kept(index)`, `constant(T)`,
  or `unavailable` entries; result-to-original identity/aggregate provenance.
  Deleted unobserved places need not be reconstructible. Names are labels.
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

## 3. Preservation is a contract, not an enum of examinations

The kernel accepts a small set of validated profiles, with extra artifact
requirements. Do not expose arbitrary combinations of permissive booleans.
Every rule has a documented applicability table and its own structural guards.

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
For duplicates, map to a representative original transition. For a simple
agglomerate store a DAG recipe such as `sequence(h, repeat(f,k))`, with original
transition IDs at leaves, rather than concatenate names. Record any pre-fired
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
be lost. This is explicit in the request, avoiding the Java coupling where the
mere presence of metadata changes the reduction policy. Later ghost support
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
for place-oriented passes, scoped by revision. A mutation either updates all
materialized affected rows or invalidates the view through one central path.
Start with invalidation and pass-level batching; incremental adjacency is a
measured optimization, not duplicated handwritten updates in every rule.

Batch independent deletions and compact with one old-to-new map. Rewrite all
matrix keys, observations, identities and metadata with that same map. IDs in
a candidate are valid only within its revision; external original IDs never
change. O(P+T) tables at preprocessing boundaries are acceptable; do not add
dense work to downstream firing loops. Tombstones and permanent stable internal
handles are optional later optimizations, not a prerequisite for this API.

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
std::optional<Match> find(const NetView<T>& net, SearchBudget& budget) const;
Estimate estimate(const NetView<T>& net, const Match& match) const; // optional
Edit<T> prepare(const NetView<T>& net, const Match& match) const;
```

These signatures express responsibilities; template placement and the exact
search cursor/result types remain design choices. No inheritance hierarchy or
runtime plugin registry is needed initially: named rule types and explicit
calls in the scheduler suffice. `find` must distinguish exhaustion from a
budget stop in its eventual result type. A match carries its source revision;
`prepare` rechecks validity if that revision changed.

The guard has two parts: eligibility for what the caller wants to preserve,
then structural premises on this candidate. Artifact maintenance is checked
before committing the edit. Both parts are mandatory even if the scheduler
already filtered the rule family. A rule may expose distinct named variants
when their proof obligations differ, rather than nested boolean switches.

A cost metric is useful but provisional. Start with a coarse search class
(local, full sparse scan, pairwise, expanding) and, where cheap to obtain,
estimated candidate work and added/removed arcs and transitions. These are
different quantities: an expensive search can produce a cheap result, and
removing a place can enlarge the transition relation. Do not collapse them
into one unexplained score. An absent estimate means unknown cost; it never
means zero. Exact measured time and actual edit deltas belong in rule stats.
The first scheduler can use fixed phases and hard budgets, then use measured
data to decide whether richer cost-based ordering pays for itself.

Each rule states: mathematical preconditions; supported profiles; observation
and artifact obligations; sparse recognition algorithm; edit; and validation
examples. Split cheap rejection tests from the semantic predicate and split
that predicate from construction. Share agglomerate construction, not an
opaque nest of `doComplex`/`doSimple` switches across different proofs.

First schedule: cheap local cleanup to stability, structural siphon/implicit
and relevance passes, non-expanding agglomeration, then budgeted broader
fusion/composition. Return to cleanup after a successful expensive rewrite.
Use actual committed change, including guard or marking changes, as progress;
net size alone is not sufficient. End at stability only if every enabled
eligible search completed without a match; report exhausted candidate budgets
as limited search, not global irreducibility.

Keep a named reference-inspired schedule for experiments, with thresholds in
one options structure. Candidate order should be deterministic (ID tie breaks)
and statistics should make schedule comparisons possible. Budget new arcs,
new transitions, candidate pairs, and total work, not only net size. Work
limits bound repeatable tests; a deadline bounds practical latency. Check
cancellation between candidate batches and before expensive allocation.

For example, post-agglomeration through initially empty p with unit arcs,
one producer h and one consumer f, h's only output p and f's only input p,
and both invisible, replaces h,f by h;f. Outside p its pre is pre(h), post is
post(f). Preconditions exclude h=f. The intermediate marking is hidden, so
this is not a next-step preservation rule. Broader weighted composition uses
k=post(h,p)/pre(f,p) only after divisibility and continuation conditions are
proved; build pre(h)+k*pre(f) and post(h)+k*post(f), omitting p, using checked
arithmetic. This formula is not a generic composition of arbitrary transitions.

## 7. Implementation stages after design review

1. Introduce request/result, revision and remapping machinery; implement
   duplicate transitions, proven constants and dead transitions. Preserve
   observations and basic counting records, with explicit refusal of any
   unsupported required artifact. Integrate one preprocessing point in
   PetriSpot, then vendor through libHSC's existing upstream-copy workflow.
2. Add empty siphons, local implicit-place detection and reachability relevance.
   Establish deadlock guards independently, especially source transitions.
3. Add the trivial chain above and simple pre/post agglomeration; compare
   schedule costs before expanding rule coverage. Decide separately whether
   original witness lifting warrants implementation.
4. Add free SCC/future fusion and broader weighted/partial agglomerations only
   with per-profile justification. Temporal, counting and dynamic-image support
   are separate milestones, not automatically inherited from reachability;
   dynamic-image support is excluded unless a concrete LTL use calls for it.
5. Add optional invariant/LP/SMT-derived facts through a separate adapter if
   useful. Structural preprocessing remains usable with no solver dependency.

Validation should target semantics rather than reproduce implementation.
For small bounded examples exhaust original and reduced reachable markings:
compare projected valuation sets, deadlock existence, or weighted counts as
the contract requires. Test near misses: read arcs, nonzero initial p,
indivisible weights, two competing consumers, visible intermediate states,
source transitions, duplicate chains, and overflow in composition. Test mixed
property supports, deletion followed by fusion, metadata reindexing and replay.
For CTL use formulas distinguishing branching, next steps, and deadlock.

Use ITS-Tools as a differential reference, not the only correctness oracle;
different reduced net sizes or rule counts are expected. Measure reduction time,
peak storage, arcs and downstream analysis cost on user-selected MCC models.
Individual diagnostics stay within the repository's 15-second bound. No tests
or benchmarks are needed for this documentation-only stage.

## 8. Review decisions

The proposed starting point is observed reachability plus a conservative
deadlock/local-counting subset, with general CTL restricted to strong local
rules. Java image support is excluded. Original witness reconstruction is
deferred; when traces are required, use only rules with an implemented lift
or run the original net. Each rule is its own file/type, with explicit guards;
cost estimates are optional scheduling advice rather than semantic premises.
The workspace is private, using ordinary dense indices per revision and batch
compaction. libHSC reuses the native API through vendoring; PNET stays an adapter.
These choices are the first implementation boundary to validate.
