# Rerun of 2026-09-06 (evening), product 202609061624 then the master walker

Submitted at 18:51 from the cluster head (`submit-2026-09-06c.sh`: QLA, SMA,
UBA, RD, QL, SM, OS, 1800 s, 4 cores, `tall`), complete at 03:00 on
2026-09-07. Logs in `/data/ythierry/MCC26run/2026-09-06c/`. The `petri64`
inside the product was replaced at 21:25 by the CI-built master `6a5d288`
(the `==` goals quested, restarts on hopeless markings, components as a
covering under a work guard, the flows boxed): QLA and most of SMA ran on the
`1624` binary, UBA onward on the fixed one. Read against
`csv/2026-09-06/` (QLA, SMA, UBA, RD of the morning) and `csv/2026-09-05/`
(QL, SM, OS).

* `total-runs.csv`: QLA 142 models better / 103 worse, answered 10.20 M ->
  9.19 M (two causes, fixed in `6a5d288`: the quest sweep where rarity won,
  and calls of 100 to 340 s lost to the components' construction, see
  HANDOFF.md 21:30); SMA 59 / 42, 4.25 M -> 3.81 M (mostly the old binary);
  UBA 178 / 56, 2.57 M -> 2.60 M (the fixed binary).
* `runs.csv`, `verdicts.csv`: RD 0 wrong, missed 9 -> 7, the three Shield
  deadlocks found, DatabaseWithMutex-PT-20 proved in 205 s; QL 1 781
  answered (1 786), missed 22 (20); SM 1 851 (1 838), missed 37 (49),
  failures 9 (14); OS 1 945 (1 942), missed 0 (1). No wrong verdict.
