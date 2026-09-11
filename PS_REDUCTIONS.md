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

## Preparation pipeline (`reduction/Pipeline.h`)

Reducing on the union support of every formula left AirplaneLD-PT-0010
untouched (89 places, 88 transitions) where ITS-Tools reached 46 places and 70
transitions: the reference substitutes the constant places into the formulas
and answers what the initial marking decides before it computes the support,
then reduces again each time a property falls. `prepare` is that loop, run by
every entry point, `--reduce` or not: constants of the current net (unchanged
by every transition, or in the greatest empty siphon) into the formulas; the
initial marking (`expr/InitialState.h`: decided formulas answered, an until
whose left side fails initially replaced by its right side, `EF p` / `AG p`
requalified as reachability kinds); decided properties reported and dropped;
reduce for the kinds and supports that remain; again while something changed.
A net declared one-safe by its NUPN structure (`SparsePetriNet::isSafe`,
cleared by the two rules that fuse places) bounds every place by one, and an
atom every value of its form decides the same way folds to a constant.

AirplaneLD-PT-0010, 15 s, verdicts all agreeing with the oracle:

| examination | answered before any engine | reduction | answered in 15 s |
|---|---|---|---|
| RC | 14 / 16 (ITS-Tools: 12) | 89 -> 17 places, 88 -> 25 transitions | 16 / 16 |
| CTLC | 5 / 16 | 89 -> 36, 88 -> 46 | 16 / 16 |
| RF | | 89 -> 54, 88 -> 88 | 11 / 16 |
| UB | | 89 -> 26, 88 -> 35 | 5 / 16 |

The reduction's own deadline is `--reductionMs`, 15 s by default, spent inside
the run's budget: ErlangenMainframeV1-PT-bP09C09 RC spends 13.7 s reducing
59 403 transitions to 52 683 and leaves the walk nothing at a 15 s cap. A
budget shared with the engines is still to design.

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
* `reach-deadlock-full.jsonl`: complete-inventory RC/RF/UB/RD campaign finished
  across all 1,681 P/T models using the fixed petri64-complete binary. There are
  zero original/reduced conflicts, 14 original and 137 reduced invocation
  timeouts, 151 model-deadline records and 46 reduction-limit reports. Stored
  comparisons record 10,050 oracle matches and zero wrong answers, but 30,112
  answers are unverified and the oracle changed during this campaign. Re-score
  against the rebuilt oracle before drawing correctness conclusions.
  BlocksWorld preprocessing reaches the allowance: a performance finding,
  not evidence of graph-rule linearity.

This bounded campaign is not a solver ranking. Unknowns are neither agreements
nor errors; absent/? oracle entries leave answers unverified. Bound lower bounds
must not exceed known maxima. Reduced runs come first within the shared allowance,
so original/reduced timeout asymmetry is not a fair performance comparison.
Most scans are sparse, but graph edge generation costs sum(|pre|*|post|), future
matching and composition searches can be quadratic, and repeated fixed points
add work. Cooperative reduction deadlines do not replace the external hard cap.

## Oracle refresh

The CI archive downloaded into the local deployment was incomplete: all ordinary
formula verdicts were `?`, and all five global-property families were missing.
This was not a withdrawal of the GPPP verdict and does not settle the disagreement.
The raw CSV and locally rebuilt archive retain TRUE for GPPP-PT-C0010N1000000000
CTLCardinality-2024-08, defended by TAPAAL and 2025GOLD; TY and PetriSpot answer
FALSE. ITS-Tools reports CC in the raw row. The counterexample investigation below resolves this particular disagreement
in favor of FALSE; the cause of the reference tools’ TRUE answers is unknown.

The complete local build in pnmcc-models-2026 produced 29,871 files, with all
ordinary families populated and all total-vector lengths preserved. Its commit
11443a6 fixes Boolean provenance parsing and stops publication on failed build
steps or missing/unfilled oracle families. CI logs supplied by the user show the
raw CSV download timing out; if retries exhausted, the old script continued
with missing input, explaining the observed artifact. The local deployment is
still the broken downloaded version; the cluster was not synchronized or deleted.
Await successful CI publication before refreshing either deployment.

The reach-deadlock-full campaign was already running during the refresh. Its
stored oracle comparisons span oracle versions (and entries unavailable during
replacement may be unverified); re-score recorded answers against one fixed
oracle before interpreting aggregate oracle-match counts. Original/reduced
pairwise comparisons remain independent of this oracle change.

## GPPP CTLC index 08: verified false

The MCC26 property is
`GPPP-PT-C0010N1000000000-CTLCardinality-2024-08`, the ninth property
(zero-based 08) of CTLCardinality, not CTLFireability. The original net with
`petri64`, seed 1 and the validation budgets answers FALSE via EXPLICIT CTL_WALK.
The ordinary root trace is blank because the initial-until failure discards
its child's evidence. Solving that AG child with the same property seed
(1 + 7919 * 8) exposes the 52-transition counterexample.

The formula has the shape `E[_3PG >= 2467342475 U AG(B)]`. Initially `_3PG=0`,
so it requires AG(B) at the initial state. In B, the outer disjunction is
`Ru5P >= Pyr` or a conjunction whose first factor is
`E4P >= 228211673 OR A[FBP >= b2 U GAP >= 2707223165]`.
The counterexample reaches Ru5P=0, Pyr=1, E4P=1, FBP=0, b2=27, GAP=1.
Both disjuncts fail: the inner until has neither its target nor its left
condition at that state. Thus B is false at a reachable state, AG(B) is false
initially, and the complete property is FALSE.

The C++ probe built with undefined-behavior sanitization reported no arithmetic
fault. Independent PNML replay in Python arbitrary-precision integers verified
that all 52 transitions are enabled and the endpoint satisfies the inequalities
above. Its maximum place marking is 4,000,000,000. This is a concrete refutation
of the consensus TRUE, not an inference from absence of a sanitizer report.
VerifyPN 4.3.0 reproduces the wrong TRUE in 0.1 s on the property alone, and
reports a token total of 9 000 000 380; Marcie rejects at parse time the 14 of
16 formulas of this instance whose constants exceed 2^31 - 1, formula 08 among
them. The oracle is patched in pnmcc-models-2026 (`install_inputs.sh`, commit
edf0e56: FALSE, defended by TY, verified from the trace), the deployed and
cluster copies with it; the trace, a standalone replay and the one-property
formula file went to TAPAAL's authors. The one-shot probes are removed.
