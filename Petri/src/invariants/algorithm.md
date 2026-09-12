# Flow elimination and optional inequality harvesting

The existing PIPE solver first computes a rational kernel basis by sparse
elimination, then optionally computes nonnegative semiflows. This document
specifies the additional, deliberately incomplete phase-1 harvesting mode.

## Contract

For the variable-by-constraint matrix M passed to InvariantMiddle, return
nonzero nonnegative integer vectors b with M^T b <= 0 (decreasing) or
M^T b >= 0 (increasing). At least one effect is nonzero. Equalities remain
in the existing basis. Each vector is GCD-normalised and deduplicated within
its direction. These are certificates, not a generating basis of the cone.
No initial marking or constant belongs to this matrix-level result.

## Isolation from the existing solver

The existing pair-returning API and its callers remain available unchanged.
A separate opt-in API returns basis, permutations, and inequalities. A compile-
time boolean selects harvesting in phase 1; the default instantiation has no
collector calls or runtime enable checks in its elimination loops. Dispatch
occurs once in the CLI. Pivot selection, heuristics, arithmetic on active
columns, culling, and phase 2 are unchanged in both instantiations.

The collector is read-only with respect to active columns. It is called only
immediately before the two existing destructive exits in eliminateRowWithPivot:
the one-sign-row batch discard and the ordinary pivot discard. There is no
hook in sparse arithmetic, row-sign maintenance, clearColumn, or phase 2.
There are no parked columns, new pivots, pair generation, or extra elimination.

## Certificate extraction

At phase-1 entry B is identity and C = A B, where A is the preprocessed
constraint-by-variable matrix. Every subsequent column operation and joint
GCD division preserves C = A B.

For a discarded column, inspect b = B[j] and c = C[j]. Ignore zero or
mixed-sign b, and zero or mixed-sign c. Orient b nonnegative, reversing c's
direction if necessary. Copy and normalise b, and insert it into that
direction's set. An unrepresentable orientation is skipped, without changing
or aborting the ordinary computation. Only accepted certificates are copied.
The scan cost is linear in the discarded pair's nonzeros; storage is linear
in the collected vectors. Duplicate columns culled before elimination are
not harvested in this first implementation.

## Coordinates and result lifetime

Preprocessing normalises constraints by positive GCD, removes zero constraints,
and merges equal positive multiples. This preserves both inequality directions.
B's row indices remain original variable indices even when candidate columns
are compacted. Inequalities bypass phase-2 permutation compression and require
no lifting. C's indices and magnitudes are preprocessed effects, so C is not
exported as an original-coordinate certificate.

The opt-in API uses the existing cooperative deadline, not the detached-thread
timeout wrapper. Certificates already collected survive an empty equality
basis or an expired deadline. Overflow in the main solver retains its existing
failure policy; no complete equality basis is claimed after such a failure.

## CLI and KERS

--collectInequalities enables harvesting alongside a flow/semiflow request.
--decreasingKERS and --increasingKERS imply collection and write ordinary KERS
matrices (original variable count by certificate count). Their option names
define the relation; the binary format and --basisKERS are unchanged. A single
export invocation must select one P/T orientation, and output paths must be
distinct. For a net, P output uses place names and b^T m0 as the constant;
T output uses transition names and states the effect direction. Pure KERS text
output describes M^T b's direction, not a fabricated constant.

## Verification

The noninterference argument is structural: both instantiations execute the
same active column operations in the same order. The observer takes const
references, copies only its accepted candidate, and writes only its own result
and duplicate sets. It never calls the sparse sum/product routines or touches
their reusable buffers. Its orientation overflow is contained locally. Phase 2
receives exactly the same B and has no knowledge of the collector.

The only enabled-mode additions are sign scans, candidate copies, normalisation,
hashing, and storage. Memory exhaustion remains possible, as for any extra
output. Disabled mode constructs no collector and instantiates no observer
calls. The existing return signatures and timeout wrapper remain unchanged.

A small internal matrix check compares equality results and validates the
certificate signs against original coordinates. It exercises a sieve, mixed
effects, positive scaling, and P/T orientation with the existing options.
Keep each invocation bounded to 15 seconds. No separate test tool or benchmark
campaign is part of this proof of concept.
