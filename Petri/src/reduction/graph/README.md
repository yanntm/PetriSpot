# Structural graph analyses

`Graph.h` supplies iterative SCC decomposition and predecessor closure on a
place dependency graph. `Dependency.h` implements StructuralReduction.buildGraph:
for each transition, connect each input place to each output place. Reachability
omits unchanged output coordinates and self edges; deadlock and the three
stutter-preserving temporal modes retain both. Duplicate edges are removed.

`Stabilizing.h` follows computeStabilizing. If every transition has nonpositive
net token effect, initially classify strictly negative transitions as stabilizing.
Then reach a fixed point: a place with no positive effect from an unclassified
transition stabilizes; any transition with a negative effect on such a place
fires only finitely often. The classification concerns eventual behavior, not
constant initial markings or transitions that can immediately be deleted.

Separate rules use these services:

* FreeSCC fuses nontrivial SCCs of unobserved unit transfer transitions, summing
  markings and incident weights. It preserves transfer transitions as self loops.
  LTL excludes it; STATESPACE requires counting transport, not yet supplied.
* PrefixOfInterest uses cyclic SCCs for DEADLOCK/SI_LTL/LI_LTL/SI_CTL and property
  support plus presets of visible transitions for REACHABILITY/temporal goals.
  In deadlock mode stabilizing transitions are excluded only while finding cyclic
  seeds; closure uses the complete graph. Temporal modes additionally include
  all inputs of transitions consuming from the prefix, then close again.
  No cyclic seed in deadlock mode proves a reachable deadlock. A source transition
  is checked first: its unconditional enabling already proves deadlock freedom.
  This makes explicit the source-transition deduction also present in Java.
* LoopBack reuses the reachability prefix computation with one inverse transition
  omitted; its rule body remains separate from ordinary prefix pruning.

Pruning reproduces dropPlaces(andOutputs=true): remove unwanted places, then
remove their former consumers only if their remaining preset is empty. The
repeated flowPT test in Java is interpreted literally, not changed to a postset
check. Source transitions with no removed input are retained.

Traversal uses explicit stacks, so long model paths do not consume the C++ call
stack. Graph construction has Java's sum(|pre(t)|*|post(t)|) edge-generation cost;
it is not universally linear in net arcs. Deadline cancellation leaves the net
unmodified during analysis.
