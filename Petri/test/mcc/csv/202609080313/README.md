# Campaign 202609080313 — CTLC, CTLF, L

The first CTL examinations at scale, and the first campaign on the
**native image** (AVX2, hence `tall%`). The ITS-Tools product `202609080313`
carries both fixes of that night: `--no-V` back in the LTSmin runners
(`7f4113e0`) and the coloured `Sort[]` registered in the closed world
(`fcdae5b8`). 1954 instances per examination, 1800 s, 4 cores. Logs in
`/data/ythierry/MCC26logs/itstools/202609080313/`.

Of the 1954 logs per examination, 58 come from the later build
`202609081701`: those instances died on the native image's closed world in
the first attempt and were resubmitted at the campaign's own budget.

|  | answered | ok | wrong | missed | none | bonus | at wall | wall h | walker h |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| CTLCardinality | 26069 | 25051 | **1** | 1731 | 3464 | 1017 | 946 | 578.4 | 44.7 |
| CTLFireability | 23682 | 22579 | **2** | 2301 | 5281 | 1101 | 1071 | 636.1 | 61.0 |
| Liveness | 1814 | 1811 | 0 | 27 | 113 | 3 | 117 | 88.1 | 25.2 |

`missed` is a value the consensus has and we do not, `none` a value nobody
has, `bonus` a value only we have. 2121 bonus values over the three
examinations is what the checker adds to the field; 4032 missed is what it
still owes, and 1505 of CTLCardinality's 1731 are simply runs at the wall.

## The three wrong verdicts

| formula | oracle | ours | run s |
| --- | --- | --- | ---: |
| `ShieldIIPs-PT-002A-CTLCardinality-2024-07` | FALSE | TRUE | 307 |
| `ClientsAndServers-PT-N0020P1-CTLFireability-2024-09` | FALSE | TRUE | 1801 |
| `FileSystem-COL-N02I10B10-CTLFireability-2024-05` | TRUE | FALSE | 767 |

None is backed by another tool, so none is an oracle to dispute. Two are
TRUE where the property is FALSE and one the converse, so no single
direction of unsoundness covers them: they are three separate hunts, and
`FileSystem-COL` is coloured, which puts the skeleton approximation on the
list of suspects for that one alone. 3 wrong in 49 565 answers.
