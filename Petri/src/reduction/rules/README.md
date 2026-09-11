# Reduction rules

Each header owns one named rule, its preservation guard, recognition and edits.
The reachability/deadlock schedule mirrors the active Java structural inventory.
SI-specific completion and validation remain deferred.

* `ConstantPlace.h`: equal pre/post rows, or initially zero with no positive
  effect; remove impossible consumers and erase redundant guards.
* `EmptySiphon.h`: greatest initially empty siphon, sparse counter propagation.
* `DuplicateTransition.h`: hash pre/post pairs; retain a named representative.
* `DuplicatePlace.h`: identical rows; the smaller initial marking controls
  enabling. Only unobserved redundant places disappear.
* `NoEffect.h`: discard identity steps for observed reachability only.
* `SinkPlace.h`: discard unobserved places with no consumers, except when STATESPACE retains coordinates.
* `PrefixOfInterest.h`: Java property/cyclic seeds and predecessor closure.
* `FreeSCC.h`: unobserved unit-transfer SCC fusion.
* `LoopBack.h`: inverse-transition exclusion followed by a fresh prefix analysis.
* `TrivialPost.h`: invisible unit one-to-one continuation, redirect and retire.
* `PostAgglo.h`: simple/general/complex F-continuation with Java weights,
  visibility and marked-place guards.
* `PreAgglo.h`: simple/complex unit-arc pre-agglomeration with invisible
  producers, divergence freedom and no competing consumers of their inputs.
* `ImplicitForkJoin.h`: a unit fork/join coordinate whose guard is entailed
  by a causal predecessor chain through the other join input (depth limit 5).

Local rules preserve next-step behavior over retained observations; constant
and siphon removal skip LIVENESS until deductions can retain dead-transition
obligations, and STATESPACE keeps token coordinates and duplicate transitions. The
reachability/deadlock scope is the validation priority. Image-dependent branches
are unavailable until image transport is implemented.

Agglomeration parity increment: pre/post recognition follows the Java simple,
non-complex and complex phases, with the same visibility, unit/divisibility,
quasi-persistence and marked-place guards (image retention remains unsupported).
`../Composition.h` prepares and deduplicates weighted transition products before
mutation, appends stable slots, and updates both sparse transposes. Post retains
the cross-product limit of 32 and complex-pass throttle of 101 applications.
Free SCC fusion uses iterative graph traversal and sums component coordinates.

Reachability/deadlock completion inventory:

* `SinkTransition.h`: invisible transitions with empty postset, reachability only.
* `ScalarTransition.h`: exact positive integer multiples of both pre/post columns;
  remove the larger transition, except for liveness. Unit nets bypass pair search.
* `RedundantComposition.h`: within equal-effect buckets keep weakest presets;
  then remove transitions dominated by a guaranteed two-step sequence. Order by
  descending arc count and stop above 20,000 active transitions as in Java.
* `BoundsDominance.h`: Java's separate positive-effect dominance method, exposed
  but not scheduled: its UpperBoundsSolver call is commented out in the reference.
* `FreeAgglo.h`: empty unobserved p, invisible unit feeders producing only p,
  unit consumers, no overlap. Simple phase also limits feeder count/input degree;
  complex phase removes these two effort restrictions.
* `PartialFreeAgglo.h`: replace eligible individual feeders by compositions;
  retain p and its consumer, require at most one unit consumer and no overlap.
* `PartialPostAgglo.h`: all consumers controlled solely by one token from p,
  with both visible and invisible choices. Compose only invisible consumers;
  retain feeders, p and visible consumers. Reachability scheduler only.
* `FutureEquivalent.h`: equal outgoing degree, greedy bijection of transitions
  modulo the two place coordinates; redirect productions and sum markings, then
  delete the redundant place and its consumers. Skip buckets >=10,000 places.
* `InitialTokenMove.h`: at final stability, pre-fire an invisible continuation
  controlled solely by an initially marked source place with one consumer.
  Discard the unusable remainder and move all complete firings in one checked
  sparse arithmetic operation, preserving Java's resulting initial marking.

These are separate rule types even where another rule may subsume them. Partial
agglomeration retains its boundary place and opposite transitions; it does not
call the full-agglomeration deletion routine. Simple and complex phases of the
same Java rule share recognition, while trivial post has its own implementation.
