# Reduction rules

Each header owns one named rule, its preservation guard, recognition and edits.
The first implementation covers local cleanup, empty siphons, reachability
relevance and restricted agglomeration. The coordinator reports the implemented
subset; the full ITS-Tools schedule remains a target, not a completed claim.

* `ConstantPlace.h`: equal pre/post rows, or initially zero with no positive
  effect; remove impossible consumers and erase redundant guards.
* `EmptySiphon.h`: greatest initially empty siphon, sparse counter propagation.
* `DuplicateTransition.h`: hash pre/post pairs; retain a named representative.
* `DuplicatePlace.h`: identical rows; the smaller initial marking controls
  enabling. Only unobserved redundant places disappear.
* `NoEffect.h`: discard identity steps for observed reachability only.
* `SinkPlace.h`: discard unobserved places with no consumers for reachability.
* `ReachabilityRelevance.h`: backward closure of observed-place effects and
  their enabling dependencies; remove objects outside that closure.
* `TrivialPost.h`: invisible unit one-to-one continuation, redirect and retire.
* `PostAgglo.h`: initially empty place, single invisible consumer controlled
  solely by that place, integral feeder weights; restricted post-agglomeration.
* `PreAgglo.h`: single invisible producer with no competing consumers of its
  inputs; unit arcs around the eliminated initially empty place.
* `ImplicitForkJoin.h`: a unit fork/join coordinate whose guard is entailed
  by a causal predecessor chain through the other join input (depth limit 5).

Local rules preserve next-step behavior over retained observations; constant
and siphon removal skip LIVENESS until deductions can retain dead-transition
obligations, and STATESPACE keeps token coordinates and duplicate transitions. The
agglomeration rules initially run only for REACHABILITY and DEADLOCK. This
restriction denotes implementation coverage, not a restriction on the Java rule.
