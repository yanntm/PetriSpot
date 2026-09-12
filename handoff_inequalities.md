# Handoff — structural inequalities and projection consumers

Current state and next actions only. Design, implementation and observations:
[INEQUALITIES.md](INEQUALITIES.md). Collector contract:
[invariants/algorithm.md](Petri/src/invariants/algorithm.md). Projection policy:
libHSC `include/hsc/linear/algorithm.md`, especially sections 4 and 6.

## Next

1. Experiments are deferred: the cluster is busy. Do not submit, deploy over
   running campaigns, or inspect the queue on this thread's behalf now.
2. When the user resumes experiments, compare default approximation with
   `--no-approx-inequalities`: newly covered places, tighter bounds, answers,
   diagram size, time and memory. Start from useful bounds queries; no broad
   internal regression campaign is needed.
3. Evaluate selective retention of bounded places before assuming that every
   extra bound should enlarge every projection. The design is documented;
   implementation still retains all certified bounded places. Keep bounds
   and the retained-coordinate mask distinct in any future refinement.
4. Let consumer gaps motivate active inequality search. Parking/postprocessing
   and kernel-assisted repairs are proposals, not implemented. Target a small
   covering family with at most n inequalities and an explicit work budget,
   not a complete cone. Preserve the cheap flow path and its pivot choices.
5. Later experiment: reconstruct safe NUPN unit bounds without supplying those
   facts to the recovery run. Distinguish direct monotonicity from bounds
   implied by larger conservation sums. No NUPN recovery code exists.
6. Lower priority: accumulating T-rythms as candidate repeatable Parikh
   segments for unboundedness. A reachable executable segment is required;
   a one-sign matrix effect alone is not an initial-marking proof.

## Operating boundary

PetriSpot's standalone collector is opt-in (`--collectInequalities`); libHSC's
approximation/deadness consumer is default-on with `--no-approx-inequalities`
as its opt-out. The existing pair-returning flow API remains available.
Neutral effects stay in equality results. Keep the collector read-only and
compile-time disabled on the legacy path; no changes to pivot heuristics.

Small examples and query inputs: `Petri/examples/inequalities/`. Their README
records the useful results and missed opportunities. The thread spans both
repositories; vendored invariant changes originate in PetriSpot.
