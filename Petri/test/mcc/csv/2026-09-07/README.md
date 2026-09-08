# Campaign 2026-09-07 — RD, QLA, LTLC, LTLF

The CI product of ITS-Tools bae71d95 (the Effort contract, spotutil, Spot 2.16)
with the petri64 of PetriSpot 3f29014, 1800 s, 4 cores, `tall%`, the Eclipse
launcher. 1953 instances per classic examination, 1681 for the total one.
Logs in `/data/ythierry/MCC26archive/2026-09-07/`, `warmups/` beside them.

|  | instances | wrong | missed | bonus | at wall | wall h | walker h |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| LTLCardinality | 1953 | 1 | 1338 | 17 | 563 | 377 | 6.3 |
| LTLFireability | 1953 | 0 | 2077 | 42 | 660 | 442 | 7.1 |
| ReachabilityDeadlock | 1953 | 0 | 5 | 9 | 42 | 37 | 5.9 |
| QuasiLivenessAll | 1681 | | | | 235 | 173 | 144.9 |

No missed value that only the contest ITS-Tools produced, on any of the three
classic examinations: the Effort glean costs nothing against the field.

QuasiLivenessAll answers 9 982 377 atoms of 19 876 189, but the mean hides the
shape: 1496 of the 1681 runs are above 0.9 completion and 62 below 0.1. Among
the runs that reach the wall having closed three quarters, that point came at
a median of 102 s -- a knee, then a residue the remaining time does not close.

## The one wrong verdict

`StigmergyCommit-PT-02b-LTLCardinality-03`, `!(G((F(p0)||X(X(F(p1))))))`:
TRUE where the consensus has FALSE, answered in 258 ms by LTSmin under
`PARTIAL_ORDER EXPLICIT LTSMIN SAT_SMT`. It is a regression -- the same formula,
net (43/685 places, 107/797 transitions) and technique tag answered FALSE on
2026-09-06. The visible difference between the two logs is the knowledge step:
on 09-06 four factoids reduced the automaton from 3 states and 4 edges to 2
states and 3 edges, and a later round found six; on 09-07 the same four factoids
reduce nothing and the later round still finds four. Weaker knowledge should
cost an answer, never produce a wrong one, so the soundness bug is in that path
and the 09-06 configuration masked it. The knowledge budget is what the Effort
contract gates: a suspect, not a proof.
