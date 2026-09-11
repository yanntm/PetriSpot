# Native structural reductions

The reachability/deadlock structural rule inventory of ITS-Tools
`StructuralReduction.reduce()` is implemented under `Petri/src/reduction/`.
This is an implementation-coverage statement, not a claim of demonstrated
performance parity or complete ITS-Tools preprocessing parity. The outer Java
solver also invokes SMT and solver-specific transformations.

## Implemented rules

Each rule has its own header, guards and documented transformation. Trivial post
agglomeration remains an independent in-place column rewrite. Simple/complex
variants of a given general Java rule share that rule's recognition code.

| Family | Native files / coverage |
|---|---|
| Place cleanup | ConstantPlace, EmptySiphon, SinkPlace, DuplicatePlace |
| Transition cleanup | NoEffect, SinkTransition, DuplicateTransition, ScalarTransition |
| Deadlock deduction | SourceTransition; cyclic-SCC absence in PrefixOfInterest |
| Graph reductions | FreeSCC, PrefixOfInterest, LoopBack; separate stabilizing analysis |
| Structural implicit places | ImplicitForkJoin, Java depth-5 causal test |
| Agglomeration | TrivialPost, PreAgglo simple/complex, PostAgglo simple/general/complex |
| Free/partial agglomeration | FreeAgglo simple/complex, PartialFreeAgglo, PartialPostAgglo |
| Future fusion | FutureEquivalent, outgoing-degree buckets and greedy matching |
| Transition dominance/composition | RedundantComposition, equal effects then guaranteed two-step sequences |
| Initial marking rewrite | InitialTokenMove, final-stability invisible continuation |
| Specialized bounds method | BoundsDominance, deliberately unscheduled as in Java's caller |

Bounds follow Java's reachability reduction configuration with objective places
protected. The separate positive-effect dominance method exists, but its sole
UpperBoundsSolver call is commented out in the reference. Activating it would
change reference behavior and requires assessing the omitted negative effects.
The upper-bound solver's external inference/SMT loop is not part of reduce().

## Scheduling and engineering

The coordinator follows the Java reachability/deadlock phase order: initial
free SCC/prefix; inner cleanup, prefix, implicit, conditional SCC, trivial or
simple post to stability; then pre, future fusion, general/complex post, complex
pre, SCC/prefix fallbacks. Siphon cleanup always follows. At stability come
redundant composition, free/partial agglomeration, then token movement. Four
consecutive transition-growing outer rounds stop expansion.

Preserved limits include implicit depth 5, redundant-composition cutoff 20,000
active transitions, future bucket cutoff 10,000 places, post cross-product
cutoff 32 when both sides branch, and 101 complex-post applications per pass.
Long transition names trigger Java's whole-net renaming after the pre phase.
Optional explicit arc/pass/time limits remain configuration controls; the arc
limit defaults to unlimited rather than silently rejecting reference cases.

The workspace maintains pre/post matrices and sparse transposes with stable
slots. Retired columns are inactive, never mistaken for source transitions.
General composition prepares and deduplicates products before clearing old
columns and appending new ones. Partial rules retain their boundary place and
opposite-side transitions. Publication compacts through existing net builders.
Core, solver implementations and libHSC are unchanged.

Details worth retaining in future audits:

* Java's current quasi-persistence helper is precisely the no-competing-consumer
  test. No stronger alternative is substituted for it.
* Prefix pruning removes former consumers of deleted places only when their
  remaining preset is empty. The repeated flowPT test in Java is preserved
  literally; it is not replaced by a postset test.
* The source-transition deduction is checked before graph analysis as well as
  during transition cleanup. Sources produce no place-graph input edges and
  must not be mistaken for evidence of inevitable deadlock.
* Stabilizing transitions are excluded only while finding deadlock cyclic seeds;
  the complete graph is restored for predecessor closure.
* Future matching retains Java's greedy matching and equalUptoPerm merge-walk
  conditions. Its unusual early exits remain an audit topic, not an invitation
  to broaden matching silently.
* Scalar comparison uses immutable pair orientation. Java's swapping of outer
  loop variables is not reproduced as mutation of the next comparison.
* Pre-firing uses checked multiplication and sparse marking updates instead of
  iterating once per token. It preserves the complete-firing result and Java's
  discarded unusable remainder.

SCC representatives and unordered candidate iteration need not yield identical
names or final net sizes. Checked arithmetic reports overflow instead of Java
integer wraparound. These implementation differences are not stronger rule
premises. No C++ speedup is claimed without measurements.

## Properties and the standalone use case

`petri64 reduce -i MODEL --props FORMULAS --output DIRECTORY` writes a matched
`model.pnet`, `properties.sexpr`, and `names.sexpr`, without solving. Supported
formulas are remapped and constant place values are substituted before ordinary
expression simplification. Boolean constants remain in the exported query file.
Constant bounds retain their form and an exact known upper bound because the
existing bound syntax has no affine constant term.

In normal MCC analysis, resolved properties are reported and dropped before
walk/CTL/LP construction; only remaining goals are explored. This connects
structural results to the existing portfolio lifecycle. A transition-free net
needs no special solver path: its coordinates are constant and formulas simplify.
Observed constant coordinates currently remain in the published net; repeating
support reduction after solved-property removal can compact the pair further.

Executable original traces and transition hints retain the original model until
lifting exists. Input PNET counting records are explicitly omitted by standalone
export. Trace hooks compile out with NoTrace. Enabled hooks currently observe
rule passes; bounded per-application captures and animated PDF output remain
unimplemented, as do image transport and full SI-specific completion.

## Validation

Validation uses existing MCC P/T archives and their supplied formulas, with
existing contest oracles, through the real CLI. Each model has a hard 15-second
shared allowance across extraction and its original/reduced subprocesses. No
cluster jobs run, no new generated tests were added, and archives are extracted
only temporarily inside the existing corpus. Results are in
`/data/ythierry/MCC26logs/local/native-reduction/`.

* All three standard binaries build. Their actual current types are int/long/
  long long; binary names do not imply a native 128-bit marking type. Existing
  sibling compiler warnings remain.
* Earlier subset: `full.jsonl`, all 1,681 P/T models, 19,714 executions,
  63,983 oracle matches, zero original/reduced verdict conflicts. One CTL oracle
  disagreement occurs in both modes on GPPP-PT-C0010N1000000000, property
  CTLCardinality-2024-08 (FALSE versus oracle TRUE), with unchanged net sizes.
  There are 35 original and 55 reduced invocation timeouts plus 94 model-deadline
  records. These are results for the earlier subset, not the completed inventory.
* Graph/general-agglomeration increment: `graph-pilot.jsonl`, 25 models,
  300 executions, 882 oracle matches, no disagreements/conflicts/errors/timeouts.
* `graph-cops.jsonl`: five CopsAndRobbers models, 60 executions, 340 oracle
  matches, no disagreements/conflicts/errors/timeouts; prefix pruning exercised.
* `standalone-db.jsonl`: DBSingleClientW-PT-d1m07, RC/RF/UB/RD exported-pair
  analysis, eight executions, 63 oracle matches, no errors or timeouts. All
  16 bound formulas resolve to the oracle value 0. The reduced UB net has zero
  transitions; the separate analysis returns in milliseconds without walking.
* `reach-deadlock-full.jsonl`: new whole-corpus RC/RF/UB/RD campaign running on
  the complete structural inventory. Its fixed binary is petri64-complete.
  The first 176 processed models have 4,164 oracle matches and no wrong
  answers/conflicts/errors, but 12 reduced and three original invocation timeouts
  plus 15 model-deadline records. BlocksWorld
  instances reach the shared allowance during preprocessing; these are concrete
  performance findings to investigate, not evidence of graph-rule linearity.

This bounded campaign is not a solver ranking. Unknowns are neither agreements
nor errors; absent/? oracle entries leave answers unverified. Bound lower bounds
must not exceed known maxima. Reduced runs come first within the shared allowance,
so original/reduced timeout asymmetry is not a fair performance comparison.
Most scans are sparse, but graph edge generation costs sum(|pre|*|post|), future
matching and composition searches can be quadratic, and repeated fixed points
add work. Cooperative reduction deadlines do not replace the external hard cap.

## Oracle refresh

The deployed oracle was refreshed from CI by removing its cached directory and
archive and rerunning install_oracle.sh. The new GPPP-PT-C0010N1000000000 CTLC
file marks all 16 properties unknown (`?`), including property 08 previously
recorded TRUE. The old campaign disagreement therefore does not establish a
regression against the current oracle. Neither oracle correctness nor absence
of overflow in PetriSpot has been proved by this refresh. The deployed format
contains no per-tool defenders for this entry.

The reach-deadlock-full campaign was already running during the refresh. Its
stored oracle comparisons span oracle versions (and entries unavailable during
replacement may be unverified); re-score recorded answers against one fixed
oracle before interpreting aggregate oracle-match counts. Original/reduced
pairwise comparisons remain independent of this oracle change.
