# Native structural reductions

Initial implementation under `Petri/src/reduction/`, following its design.
The module is owned by PetriSpot; reference repositories are read-only. No
cluster execution is used. Local model validation has a 15-second hard limit.

## Scope

First increment: isolated workspace, named goals, coordinator, 11 rules covering
local cleanup, empty siphons, fork/join implicit places, reachability relevance,
trivial post and restricted pre/post-agglomeration. All
ITS-Tools goal names remain recognizable; complete reference rule/schedule
parity is not claimed. High-debug rendering, image support, counting metadata
transport and executable trace lifting remain separate increments.

## Engineering questions

The existing `SparsePetriNet` has no bulk replacement preserving names, so the
module publishes through its existing named builders. This adds one sparse
linear pass without changing core APIs. Retired transition slots are explicitly
inactive; clearing their columns never creates source transitions.

Observed constants initially remain as isolated places. That preserves the
existing property coefficient width without needing wide constant substitution.
Original executable witnesses and Parikh hints need lifting across rewriting;
the CLI must retain the original net when these artifacts are requested.

`--reduce` is opt-in and runs after property parsing but before existing walk,
CTL or LP compilation. Deadlock queries use their own reduction goal. Mixed
deadlock/state queries and general CTL use the conservative strong local subset.
The original input is kept for exports/invariants and for trace/hint requests.
Only CLI Options, WalkDriver and LpDriver require code edits outside the module;
core, expression trees, solver implementations and libHSC are unchanged.

## Validation

The maintained suite is `Petri/test/reduction/`: existing MCC instances and
their supplied formulas through the real CLI, original and reduced, against
the existing contest oracle. The small finite-state check from bring-up is
retained as a supplementary diagnostic; validation effort goes into MCC models
and formulas rather than expanding generated tests.
Each model is a separate worker capped at 15 seconds across archive extraction
and all of its analysis subprocesses. Archives are extracted temporarily inside
the corpus and cleaned after use. No cluster runs or benchmark copies in this
repository. Raw results live under
`/data/ythierry/MCC26logs/local/native-reduction/`.

* All three binaries (`petri32`, `petri64`, `petri128`) build. Compiler output
  contains existing sibling warnings; no reduction-header warnings were found.
  Actual build types are int/long/long long; the binary names do not imply a
  native 128-bit marking type in the current CMake configuration.
* `pilot.jsonl`: 25 P/T instances, 300 original/reduced CLI runs across RC, RF,
  UB, CTLC, CTLF and RD. 805 oracle matches, zero oracle disagreements,
  pairwise conflicts, errors or timeouts. This pilot preceded addition of the
  restricted pre and implicit fork/join rules.
* `lp-pilot.jsonl`: 25 P/T instances, 150 original/reduced LP runs over RC/RF/UB.
  173 oracle matches (85 original, 88 reduced), zero oracle disagreements,
  conflicts, errors or timeouts.
* `full.jsonl`: whole P/T corpus pass in progress; final coverage and exceptions
  will be recorded after completion. The binary stays fixed throughout it.

These are bounded validation runs, not a solver ranking or proof of all rule
implementations. Unknown answers count neither as agreement nor as disagreement.
An oracle `?` or absent entry leaves an emitted answer unverified. UB lower
bounds are checked not to exceed known oracle maxima. Reduced runs precede
original runs within each shared model allowance; timeout asymmetry therefore
precludes a fair performance/coverage comparison from this campaign alone.
Reduction summaries distinguish preprocessing time from analysis/parse time.

## Reference observations and next coverage

Existing ITS-Tools RC logs under `itstools/2026090600/RC` identify reduction-rich
examples: `OAR.1333844.stdout` (CopsAndRobbers), `OAR.1334338.stdout`
(HirschbergSinclair), `OAR.1333818.stdout` (CloudReconfiguration). Keep these
as external references when bringing additional rule variants to parity.

The coordinator currently follows cleanup-to-stability then siphon/implicit,
trivial post, single-consumer post, single-producer pre, returning to cleanup
after progress. It does not yet claim the reference's complete nested schedule.
Missing families include free SCC/future-equivalence fusion, broader and partial
agglomerations, redundant compositions/scalar multiples, token pre-firing and
the reference's temporal/deadlock SCC reasoning. Their existing applicability
conditions and model-specific performance limits need explicit transcription.

STATESPACE currently keeps every place and duplicate transition rather than
reconstructing token/arc metadata. LIVENESS skips constant/siphon removal so
dead-transition obligations are not lost. No API yet transports PNET counting
records or Java images. Pass-level trace-policy hooks compile out under
`NoTrace`; application-level local captures and PDF rendering remain pending.

Most current scans are linear in sparse structure, but repeated fixed points,
hash collision buckets and depth-bounded causal search are not a universal
linear-time guarantee. `--reductionMs` is checked cooperatively between batches;
the external per-model timeout is the hard limit. Weighted agglomerates are
prepared with checked arithmetic before mutation; an overflow currently raises
an explicit error rather than returning a partially reduced result. Improving
that recovery belongs in a later increment.
