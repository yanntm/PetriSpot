# Baseline of the reachability and LTL examinations, 2026-09-06

First campaign of RC, RF, LTLC, LTLF (`submit-2026-09-06b.sh`, from 03:32 on
2026-09-06, product 202609060003, 1800 s, 4 cores); CTLC and CTLF were
submitted and deleted unrun. Logs in `/data/ythierry/MCC26run/2026-09-06/`
(stdout) and whole in `/data/ythierry/MCC26archive/2026-09-06-baseline/`.

| examination | answered | ok | wrong | missed | bonus | failures |
| --- | --- | --- | --- | --- | --- | --- |
| ReachabilityCardinality | 30 562 | 30 427 | 2 | 381 | 133 | 21 |
| ReachabilityFireability | 30 048 | 29 976 | 1 | 728 | 71 | 20 |
| LTLCardinality | 28 932 | 28 916 | 0 | 1 457 | 16 | 13 |
| LTLFireability | 27 413 | 27 395 | 0 | 2 548 | 18 | 21 |

The three wrong values are the two known ones: NeoElection-PT-3 and PT-7 in
RC (the walker's FALSE read as a witness, fixed in ITS-Tools `3aeb7f95`) and
LastZero-COL-N20 formula 13 in RF (the skeleton's enabling polarity, fixed in
`d26cddfa`); both fixes are in product 202609061502 and later.
