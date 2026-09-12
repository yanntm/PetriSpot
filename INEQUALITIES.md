# Structural inequalities: computation and first consumers

This thread extends the existing flow computation with a small set of
structurally monotone sums. The immediate objective is useful coverage for
bounds and behavioural understanding. A complete generating set of a cone is
not the objective. The first implementation is an opt-in POC inside PetriSpot;
consumer experience should determine whether a more active search is worth it.

## Matrix contract

The public matrix API takes M with variables as rows and constraints as columns.
It collects nonzero nonnegative integer vectors b satisfying either:

* decreasing: M^T b <= 0, with at least one nonzero effect;
* increasing: M^T b >= 0, with at least one nonzero effect.

Zero effects are ordinary flows/semiflows and are not reported again as
inequalities. Mixed-sign effects are not monotone certificates. Coefficients
are GCD-normalised. Each reported result is sound; absence proves nothing.

For P vectors of C = Post - Pre, these certify b^T m <= b^T m0 or the reverse.
A decreasing vector bounds every place p with b_p > 0 by floor((b^T m0)/b_p).
An increasing sum need not bound any one place below: tokens can move between
its members. Lower bounds initially equal to zero still express a monotone
trend, although the static lower bound alone may add no information.

For T vectors, use the transposed incidence matrix: x >= 0, Cx of one sign.
We call these transition rythms. Cx >= 0 and Cx != 0 is an accumulating rythm;
Cx <= 0 and Cx != 0 is a consuming rythm. They are algebraic multisets, not
necessarily executable from the supplied initial marking. A reachable,
executable self-covering segment with a positive gain witnesses unboundedness;
the structural vector alone does not establish its reachability. For ordinary
P/T nets, an accumulating vector can be executed and repeated from some
sufficiently large initial marking.

## Implemented POC

Local commit c533c17 introduces the collector, CLI/API and five small examples.
The implementation and isolation contract are described in
[invariants/algorithm.md](Petri/src/invariants/algorithm.md).

* `Inequalities.h` owns the read-only collector and its two output matrices.
* `InvariantCalculator.h` has exactly two observation boundaries in phase 1:
  immediately before the one-sign-row batch discard, and before clearing an
  ordinary pivot. No callback is inserted in sparse arithmetic or phase 2.
* `InvariantMiddle::computePInvariantsWithInequalities` returns `basis`,
  `permutations`, and `inequalities.decreasing/increasing`. The existing
  pair-returning entry points remain available unchanged.
* `--collectInequalities` enables harvesting alongside P/T flows or semiflows.
  `--decreasingKERS` and `--increasingKERS` imply collection and export ordinary
  KERS matrices. P and T exports use separate invocations.
* P text output uses names and initial-marking constants, with a monotonicity
  annotation. T output uses transition names and the nonzero effect direction.
  The matrix API and KERS contain coefficients, not marking-dependent constants.

Phase 1 maintains C_work = A B. A discarded B column is already expressed in
original variables; its C_work column supplies the effect sign. Constraint
normalisation uses positive GCD scaling and positive-multiple deduplication,
so it preserves those signs. Candidate-column compaction does not change B's
variable indices. Inequalities bypass the later permutation compression.

The legacy solver instantiation compiles out harvesting; the collector cannot
change active matrices, pivot selection or heuristics. Enabled mode adds scans
of discarded sparse pairs and copies/hashes only qualifying candidates. It
does not park candidates, search combinations, or harvest overwritten active
columns. There is no additional flow run merely to produce the certificates.

The new API uses cooperative deadlines. Collected certificates remain useful
even when the equality basis is empty or phase 1 cannot finish. They must not
be suppressed solely because the equality result is empty. Collector-only
orientation overflow skips that candidate; main computation overflow follows
the solver's existing failure policy. Consumers must check constant arithmetic.

## Small models and what they taught us

The executable inputs and complete observed P/T results are in
[examples/inequalities](Petri/examples/inequalities/README.md). All arcs have
weight one. For net 2 we instantiate N = 5.

| Net | Behaviour | P-flow inequalities actually collected |
| --- | --- | --- |
| 1 | p0 -> p1; p1 -> p0 + p2; initially p0=1 | p2 >= 0 |
| 2 | p0 -> p1; p1 + p2 -> p0; initially p0=1, p2=5 | p2 <= 5; p0+p2 <= 6 |
| 3 | p0 -> p1; p1 -> nothing; initially p0=1 | p0 <= 1 |
| 4 | nothing -> p0; p0 -> p1; initially empty | p1 >= 0; p0+p1 >= 0 |
| 5 | nothing -> p0; p0 -> p1; p1 -> nothing; initially empty | none |

Nets 1 and 2 also have p0+p1=1. Net 5's t0+t1+t2 cycle is neutral and stays
in equality results. Flow and semiflow runs encounter different intermediate
candidates, hence can collect different inequalities despite agreeing on an
equality basis for these examples. P/T flows and semiflows were exercised;
net 2's decreasing vectors were also exported and decoded as KERS.

Net 3 misses p0+p1 <= 1. Its p1 candidate has effects (1,-1), and the sink
equation removes it before the later p0 direction, effects (-1,0), can repair
it. Keeping both would permit (1,-1)+(-1,0)=(0,-1).

On transitions, nets 1 and 2 miss t0+t1, whose effects are (0,0,1) and
(0,0,-1). Individually the transition effects have mixed signs, and the
equality-oriented order discards the useful direction too early.

Another lost opportunity is an active candidate overwritten before it reaches
a discard boundary: net 5's source t0 already has a positive effect initially.
Observing such candidates and parking discarded candidates are distinct possible
extensions. Neither is implemented or required of this POC.

## More active search: agreed objective, not implemented

Keep cheap harvesting as one level. A portfolio may explicitly trigger a
more active pass when decreasing coverage is low, or when harvested vectors
suggest promising structure. The active result should have at most n
inequalities for n variables, prioritising additional useful coverage and
simple coefficients/supports. Completeness of the cone is secondary. The
cheap POC currently has no explicit output-limit option; the active contract
must state and enforce its limit rather than inherit an accidental count.

Track decreasing and increasing coverage separately. A place's membership in
a decreasing sum supplies an upper bound; membership in an increasing sum
describes a different guarantee. A positive equality already provides bound
coverage, whereas a mixed-sign flow does not necessarily do so. Coverage alone
is not sufficient as a score: adding all decreasing certificates creates one
large covering sum, potentially losing tight individual bounds and independent
trends. Prefer certificates that add coverage, sharpen bounds, or simplify a
previous certificate. Work must also be budgeted; small output does not imply
a cheap search.

The active API would accept the original matrix, known equality results, known
inequalities, optional target variables, and effort/output limits. No flag or
specific active algorithm has been implemented yet.

### Candidate algorithms considered

1. **Parking during elimination.** Retain discarded coefficient/effect pairs
   outside active pivot selection. Later pivots also update affected parked
   candidates. Active directions are zero on previously eliminated equations,
   so retained effects there are not disturbed. Batch discards need an explicit
   extra cancellation step: they currently perform no pivot arithmetic.
2. **Postprocessing retained directions.** Save discarded pairs, then replay
   later cancellation directions afterwards. This separates the extra work
   from the ordinary computation. The final kernel alone cannot recover the
   discarded nonzero effects. Retain only as much information as the budget
   and intended coverage justify.
3. **Dedicated bounded heuristic.** Start from coordinate candidates and their
   effects; collect immediate signs, repair mixed effects through selected
   combinations. Limit continuations instead of generating all sign pairs.
4. **Targeted linear feasibility.** For a desired place p, seek x >= 0,
   Ax <= 0, x_p >= 1 (or reverse the effect sign). This targets coverage without
   cone enumeration. Zero-effect answers belong with equalities. Rational
   candidates require an exact integer certificate before publication.
5. **Complete cone computation as a reference target.** Ax <= 0, x >= 0 is
   equivalent to Ax+s=0, x,s>=0. Its graph has explicit basis [I;-A], allowing
   a nonnegativity phase without preliminary Gaussian elimination. Direct
   double description instead retains the satisfying side of each inequality
   and generates boundary combinations. Both can grow exponentially and are
   not the desired default. Support-minimality pruning in x alone is invalid;
   the effect/slack coordinates or active constraints matter.

### Reusing known equalities

Given AK=0 and a candidate v with a desirable effect, seek x=v+K lambda >= 0.
This repairs coefficients without changing Ax=Av. Lambda is unrestricted for
a linear kernel basis. It cannot repair a mixed effect; other directions are
needed for that. A full kernel basis may reduce a search to rank(A) effect
dimensions, but the coordinate change can destroy sparsity. Known semiflows
need not span the whole kernel: distinguish a full basis from a collection of
valid equalities in the active API. Kernel-assisted repair is useful even when
we choose not to form a quotient representation.

## First consumers: inspected integration points

### Later direction: transition rythms as unboundedness candidates

A nonnegative T vector x with Cx >= 0 and (Cx)_p > 0 is a candidate Parikh
vector for a repeatable segment growing place p. It can guide the search for
a firing order and a reachable marking from which that segment executes.
Once an executable segment from m to m+Cx has been found, ordinary P/T
monotonicity makes it repeatable and proves p unbounded. The matrix vector
alone is not that proof: its enabling resources or ordering may be unavailable
from the initial marking. Use the existing walks to establish the prefix or
realise the segment, rather than declaring unboundedness from a positive
effect. This could find proofs faster than unguided self-covering walks. It
is a later, non-priority direction; no finisher search is implemented.

### Experiment direction: recovering NUPN unit constraints

A useful future experiment is to hide the NUPN facts from the computation,
then compare the recovered constraints with each safe unit's
sum_{p in unit} m_p <= 1. Distinguish a direct nonincreasing indicator sum
from a bound implied by a larger decreasing sum or a conservation equality
involving places outside the unit. A safe unit can become empty and then
occupied again: its local sum is bounded by one without being monotone, so
the indicator need not satisfy C^T b <= 0. Enabling conditions can also prove
unit safety beyond this matrix abstraction. Measure units certified, places
newly bounded by one, and the size/coefficients of the explanations; do not
expect to reconstruct every NUPN fact. Avoid circular evidence: supplied safe
tags and unit constraints must not participate in the recovery run. This is
an experiment proposal only; no recovery algorithm or experiment is implemented.

### libHSC projection / invariant approximation

`tools/pn_approx.hh::approx_facts` already calls `hsc::petri::pflows` on the
full net. It derives individual bounds from positive equality vectors, then
adds safe-unit bounds and structural zeros. `build_approx` constructs a box F,
flow equalities, unit at-most constraints, and their intersection S. This is
the most direct first consumer: unlike the full state equation, that current
approximation need not already imply our monotonicity inequalities.

Enabled by default in libHSC's approximation/deadness pass. Use
`hsc-pn --no-approx-inequalities` to disable it; `--approx-inequalities`
remains accepted for explicit enabling. Standalone PetriSpot flow scenarios
remain opt-in, as before.

* The native bridge (`include/hsc/petri/invariants.hh`,
  `src/petri_invariants.cc`) has `pflows_with_inequalities`, carrying separate
  equality, decreasing and increasing lists. The existing equality API for
  shape/decomposition callers is unchanged. Required upstream files were
  copied through libHSC's vendor script, which now supports an explicit subset.
* Decreasing vectors feed the existing per-place upper-bound calculation:
  min(current bound, floor(constant/coefficient)). Do this before model/domain
  construction, as approx_facts already does for equalities.
* Decreasing sums are at-most filters on S; increasing sums use the existing
  signed at-most constructor with both coefficients and constant negated.
  All relation types stay distinct. Bounds alone lose useful correlations.
* Facts are computed before projection. Only inequalities whose whole support
  survives are kept and reindexed. Decreasing supports are bounded by their
  own certificates. Increasing supports may be removed; their terms are never
  silently dropped. Forward induction suffices for these constraints; they
  are not asserted to be invariant under reverse firing.
* Statistics report retained inequalities, additional bound coverage and
  tighter bounds. Disabled, the approximation still calls the original bridge.
* Bound queries over retained places can now maximise their sum over S in the
  inequality-enabled path (`pn_approx_bounds.hh`). This yields an upper bound. It becomes
  an exact numeric FORMULA only when attained by the initial marking; otherwise
  the statistic is retained in the output and the query stays open. No extra
  reachability or flow computation is introduced.

The user explicitly authorized this libHSC consumer. It does not change
PetriSpot's LP or walker scenarios and does not implement active search.

### Consumer observations

Using `--approx 2 --approx-only --printUnknown`, comparing harvesting enabled
and disabled (the original observations preceded making it the default):

* Net 2: coverage rises from 2/3 places to 3/3. The approximation now has 12
  markings over all three places; the old 2-marking set projected p2 away, so
  those cardinalities are not comparable as a measure of precision. All five
  resource queries formerly UNKNOWN are answered: p2<=5, p0+p2<=6, impossibility
  of p2>5, and exact maxima 5 and 6. Both maxima are attained initially.
  For p1 the approximation upper bound is 1 and the initial lower bound is 0;
  that maximum remains UNKNOWN in this pass rather than being asserted attained.
* Net 3: coverage rises from 0/2 to 1/2. Its p0 projection has two markings,
  proves p0<=1 and maximum p0=1. The total bound remains UNKNOWN because p1
  is still uncovered; the example exposes the collector's missing total sum.
* Net 1 gains no upper-bound coverage; its increasing p2 certificate cannot
  survive removal of p2. Nets 4 and 5 likewise gain no bounded projection.

The newly enabled one-place projection exposed a pre-existing boundary in
libHSC's linear-set surface reader: it reads domains from product arcs, not
a bare leaf root, producing an empty filtered set for the collapsed shape.
The inequality-enabled projection now emits `(spine place)` for that case. Its unit tail
adds no variable or state. Net 3's approximation is now the expected nonempty
two-marking set; a query containing its initial marking is not falsely refuted.
No calculus-core refactor was needed.

### Cost of making the consumer default

Harvesting reuses the same flow computation, but is not literally free: it
scans discarded pairs and stores certificates. Under a deadline this overhead
can reduce how much elimination finishes. More importantly, stronger bounds
can retain additional places and larger domains in the projection, while
additional filters cost diagram work. Precision improves for the same facts,
but total runtime and memory can improve or worsen. The opt-out keeps a direct
comparison and fallback. No active search is enabled by this default.

### PetriSpot LP and bounds

`cli/LpDriver.h::runLp` currently builds the full state equation; it does not
compute flows. Adding our inequalities as rows cannot strengthen that exact
relaxation: they already follow from m=m0+Cx, x>=0. Avoid paying for a flow
computation merely to add redundant rows.

Useful roles are instead cheap per-query bounds, exact certificates retained
when a numerical solve times out, or closing a bound without a solve when its
known upper bound equals a reachable lower bound. A matching decreasing sum
is maximised at the initial marking itself. More generally a covered target
can be upper-bounded by a suitable nonnegative multiple of a certificate,
though this may be loose. Full LP dual certificates can explain stronger
bounds and should not be confused with the small harvested family.

Any LP consumer should accept already computed results rather than silently
rerun flows per property. Bound-query consumers must preserve the distinction
between a certified upper bound and an attained maximum.

### Existing walk component construction

`cli/WalkDriver.h::buildComponents` is the only non-CLI-analysis flow call in
PetriSpot itself. It computes semiflows for component-guided walking. It is
not currently a bounds consumer. Do not append inequalities to the semiflow
list: conservation-dependent component logic has a different contract. An
independent structural-bound store could consume them later if the walker or
portfolio needs it; this is not grounds to change the current components.

## Next work

Use the projection consumer's coverage and missed opportunities to decide
whether parking/postprocessing or a targeted active solver is worth adding.
The NUPN recovery comparison remains a documented experiment direction, not
an implementation. Keep the cheap default and existing equality scenarios
unchanged; keep internal checks modest and focus on useful bounds.
