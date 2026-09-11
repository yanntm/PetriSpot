# ITS-Tools structural reduction review

Inspection of the local working files, 2026-09-11. Method names are the stable
navigation anchors; line numbers below describe this inspected snapshot.
Reference trees were read only. This is a code review, not a proof audit or a
performance measurement. Reasons labeled as interpretations are inferred from
the implementation; no external literature claims are made here.

## Source map

All Java paths below are relative to `~/git/ITStools`.

* `fr.lip6.move.gal.structural/src/fr/lip6/move/gal/structural/StructuralReduction.java`:
  about 3,360 lines, read in full. Net ownership, rule bodies, schedule, graph
  construction, images, mutation, deductions and diagnostic output coexist.
* Same folder: `SiphonComputer.java`, `NetBlocks.java`, and
  `SparsePetriNet.readFrom` (482 onward), read for the adjacent contracts.
* `pnmcc/fr.lip6.move.gal.application.pnmcc/src/fr/lip6/move/gal/application/solver/ReachabilitySolver.java`:
  `applyReductions` (688 onward) and the beginning of its SMT escalation.
* Same solver subtree: `ltl/LTLPropertySolver.buildReduced` (511 onward), for
  atomic-proposition support and `keepImage` in dynamic product approaches.

PetriSpot files inspected: `core/SparsePetriNet.h`, matrix mutation APIs,
arithmetic helpers, `expr/Property.h` and linear atom representation,
`io/PNET.md` and `PNETIO` block API, CLI property loading. libHSC's
`include/hsc/petri/README.md` establishes upstream vendoring, and
`tools/hsc-pn.cc` shows net/block loading and counting consumers. Its README's
claim that the invariant solver is not vendored conflicts with a later entry
describing a vendored calculator; no design here depends on that claim.

## Rule inventory

| Java anchor (line) | Actual responsibility | Proposed home |
|---|---|---|
| `ruleReduceTrans` (874), `ensureUnique` (1002) | No-effect/sink removal, exact duplicates, scalar multiples, source-transition deduction | Local rules; split recognition from policy |
| `ruleReducePlaces` (1088), `computeConstants` (1327) | Constants, empty places, sinks, duplicate places, pre-firing, loop-back relevance | Several local rules and relevance |
| `ruleImplicitPlace` (1384), `inducedBy` (1543) | Restricted unit-weight fork/join causal implication | Implicit |
| `ruleTrivialPostAgglo` (1574) | In-place one-to-one chain shortcut | Post agglomeration fast path |
| `rulePostAgglo` (1677) | F-continuation, weights, visibility and branching guards | Post agglomeration |
| `rulePreAgglo` (2099) | Quasi-persistence, divergence, visibility and unit arcs | Pre agglomeration |
| `ruleFreeAgglo` (796), partial variants (550, 681) | Reachability-oriented relaxed/partial composition | Free and partial agglomeration |
| `agglomerateAround` (1977) | Delete/append, weighted composition, deduplication, transpose repair | Compose + workspace |
| `ruleRedundantCompositions` (441) | Same-effect dominance, then replacement by a two-transition sequence | Dominance/composition |
| `ruleRedundantCompositionsBounds` (377) | Positive-effect dominance for bounds; outside `reduce()` | Separate specialized rule |
| `ruleFusePlaceByFuture` (2293) | Match outgoing futures modulo place permutation; sum places | Fusion |
| `ruleSymmetricChoice` (2485) | Older specialized matching, call commented out | Do not activate by transcription |
| `findSCCSuffixes` (2650), `computeSafeNodes` (2774) | Query-dependent dependency closure and deadlock reasoning | Relevance |
| `findFreeSCC` (2861) | SCC of unobserved unit transfer transitions, aggregate tokens | Fusion |
| `computeStabilizing` (3148) | Monotone token-effect reasoning to exclude infinite firings | Graph analysis |
| `dropPlaces`, `dropTransitions`, `dropPlace`, `removeAt` | Matrix/name/marking/support/index edits | Workspace |
| `createSumOfVars`, `fusePlaces` | Explicit model transformations outside main schedule | Later separate transformations |
| `abstractReads`, `applyCausalDecomposition` (3253) | Abstraction/experimental decomposition outside main schedule | Exclude from equivalence reducer |

The enum contains DEADLOCK, REACHABILITY, SI_LTL, LTL, LIVENESS, STATESPACE,
LI_LTL, SI_CTL, NONE. Individual methods contain additional mode gates: reaching
a method in the scheduler does not mean that method applies in every mode.

## The actual schedule

`NONE` exits. `STATESPACE` takes a separate short route: place cleanup without
siphons or token movement, transition cleanup, redundant compositions, free
SCC only when a counting record exists, then place/transition cleanup again.
It does not run the general fixed-point loop.

For the other modes:

1. Free SCC fusion; transition cleanup if it changed anything; SCC suffix
   pruning. SI_LTL, LI_LTL and SI_CTL additionally try place cleanup with
   token movement before the loop.
2. Inner fixed point: place cleanup without siphons/movement, transition
   cleanup, SCC suffixes, structural implicit places. If there was progress,
   retry free SCC. Try trivial post-agglomeration, and only if it did nothing,
   simple post-agglomeration. Repeat while the accumulated count is positive.
3. Try simple pre-agglomeration. If it did nothing, try future-equivalent
   place fusion. Each following step is gated on no accumulated progress:
   non-complex general post, complex post, complex pre, free SCC, SCC suffixes.
4. Always try place cleanup with empty siphons, without token movement.
   Only at stability try redundant compositions, then (reachability only)
   simple free, complex free, partial free. Then partial post for reachability
   or selected temporal modes without `keepImage`. Finally try token movement.
5. Restart the global loop if anything changed, unless more than three
   consecutive rounds increased the transition count (stop after the fourth).

This is a strategy, not a canonical reduction normal form. The likely rationale
is to expose cheap simplifications before multiplying transitions and to let
cleanup absorb the results of expensive rewrites. Exact ordering has no
standalone soundness role **provided** each rule checks all premises on the
current net. Recognition restrictions can still encode essential soundness
conditions and must not be removed as mere tuning.

There is another scheduler outside this class. `ReachabilitySolver.applyReductions`
alternates the structural routine with first-pass arc/safety-triggered dead
transition checks and optional SMT rules. It delays state-equation reasoning,
caps that case at roughly 10k places and 10k transitions, and uses
`lastReduction` to avoid cycling through unproductive stages. Porting only
`reduce()` will not reproduce the complete ITS-Tools preprocessing pipeline.

## What is heuristic, and what is a premise?

| Choice | Classification and consequence |
|---|---|
| Implicit search depth 5 | Incompleteness/effort limit; deeper proof search is possible |
| Skip redundant composition above 20,000 transitions | Cost cutoff for potentially quadratic work |
| Skip future-equivalence degree buckets of size >=10,000 | Cost cutoff, not evidence of inequivalence |
| Reject post cross-product >=32 when both sides have multiple members | Growth heuristic; distinguish from weight/choice guards |
| Stop complex post after `total > 100` | Pass throttling (can apply 101), not a theorem |
| Pre with `keepImage` rejects >=4 consumers | Growth heuristic documented with a model example |
| Stop after four consecutive transition-growing rounds | Global blow-up brake; not convergence |
| Rename when a transition name reaches 1024 characters | Representation workaround for concatenated histories |
| Sort composition candidates by descending arc degree | Search priority; interpretation: eliminate bulky composed transitions first |
| Single consumers/producers in several rules | Sometimes only cheap-first restriction; sometimes continuation, conflict or branching premise: classify per rule |
| Empty eliminated place, integral weight ratio, no producer/consumer overlap | Semantic premises for the corresponding composition |
| Invisible side, SI_CTL's single continuation, mode exclusions | Preservation premises, not optimization knobs |

The 30-second messages within post-agglomeration are progress reporting, not
timeouts. The routine has no general deadline or cancellation mechanism.

## Optimizations worth retaining

Sparse pre/post columns make transition tests merge operations on their arcs.
Transposes give place neighborhoods; several methods postpone their creation
until a candidate exists. Composition narrows second-transition candidates to
consumers of places fed by the first transition. Duplicate detection hashes
pre then post columns; effect dominance maintains weakest preconditions inside
effect buckets. Future matching first groups places by output degree.

Removal runs in descending index order to keep pending indices valid. Batch
row deletion avoids repeatedly shifting all columns. Agglomeration snapshots
old columns, creates a temporary batch, deduplicates it before insertion, and
patches already materialized transposes. Trivial post replaces a column directly
and delays cleanup rather than building a Cartesian product. `maxArcValue`
avoids scanning for scalar multiples on unit-weight nets.

These are useful algorithmic choices, not just Java idioms. Boxed collections,
stream grouping, repeated full transposes, and temporary expression trees are
implementation choices we need not reproduce. Preserve sparse locality first;
measure whether maintaining both orientations beats rebuilding at pass boundaries.
The degree filters bound some pairwise searches but do not make them linear.
Likewise the place dependency graph can expand a wide transition quadratically.

`SiphonComputer` starts with all initially empty places and repeatedly removes
outputs of a transition that feeds the candidate but consumes from none of it.
The remaining greatest initially empty siphon stays empty forever. Its repeated
scans can become a sparse worklist with counters of remaining candidate inputs;
that optimization should preserve the same fixed point and include source
transitions from the outset.

## Refactoring and correctness audit points

These are static inspection findings; no failing execution was produced in
this stage. They justify focused tests rather than claims about contest results.

* `ruleReducePlaces` returns removed-place count, although it can delete
  transitions, clear guards, or move tokens as well. `reduce()` mixes boolean
  SCC progress, deletion counts and application counts; some increments are
  overwritten or added inconsistently. Use a committed-change flag for
  termination and separate counters for reporting.
* `testModuloIsomorphism` caches `coli` before swapping `ti/tj` and `vi/vj`;
  it also mutates outer-loop locals across inner comparisons. Re-express this
  using an immutable pair and test both weight orders before reusing it.
* `dropPlaces(andOutputs)` repeats the `flowPT` emptiness test twice in its
  filter. Its comment and predicate need reconciliation; do not guess that the
  second occurrence should be `flowTP` without checking relevance semantics.
* Protected constants clear columns in transposed views, but reconstruction
  is conditional on deletion counts. This is another reason to track edits,
  not only removals, and assert view consistency.
* Several mutations bypass `dropTransitions`/`NetBlocks`; `clone()` does not
  copy blocks, and the constructor's block-copy path is specific to
  `SparsePetriNet`. Uniform metadata maintenance is not established by the
  presence of the helper alone. Native edits must have one commit path.
* `findFreeSCC` aggregates markings while other fusion paths explicitly clear
  `isSafe`; it does not visibly do so here. A component of individually safe
  places can have an aggregate marking above one. Recompute/invalidate facts.
* `image` begins as original variable expressions and some transformations
  add sums or offsets. `keepImage` also changes which consumer transitions
  survive agglomeration: it is a dynamic-product mechanism, not a logging flag
  or a complete original-trace lift. Free SCC and other paths need separate
  examination before asserting any general image invariant.
* `SparsePetriNet.readFrom` remaps properties through place names. Native
  identities should be numeric and explicit, including constants and loss of
  information; concatenated transition names are not a provenance format.
* The scheduler itself comments that partial post is “almost legitimate” for
  SI_LTL, yet invokes it for selected temporal modes. This is a proof-review
  item, not permission to ship that temporal eligibility unchanged.
* SCC deadlock deduction is attempted before `ruleReduceTrans`'s source test.
  A source transition is always enabled and prevents deadlock; a graph with
  no cyclic place SCC must not overrule that fact. Check standalone `reduce`
  on a source-only net rather than rely on caller preprocessing.
* `abstractReads` removes guards and is an abstraction. Causal decomposition
  warns on a marked place but proceeds; its per-feeder copies also need a
  separate argument for consumers requiring tokens from multiple feeders.
  Neither belongs in the default equivalence-preserving rule registry.

The useful split is therefore threefold: sound recognition with a stated
relation, centrally checked mutation with artifact maintenance, and an
independently tunable scheduler. This retains the sparse algorithms while
making both correctness and performance decisions inspectable.
