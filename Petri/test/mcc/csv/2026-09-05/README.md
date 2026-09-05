# Campaign of 2026-09-05 — global properties and UpperBounds

ITS-Tools (`BK_TOOL=itstools`, no reducer prefix) on the whole MCC 2026 set,
six examinations, one `oarsub` per pair, 1800 s per test, 4 cores, all on the
`tall` nodes of `cluster.lip6.fr`. Builds `202609051237` (the first jobs) and
`202609051349` (the rest). Logs live in `/data/ythierry/MCC26run/`, rsynced
from `~/MCC26/MCC-drivers/{L,OS,QL,RD,SM,UB}`.

1953 instances, but twelve of the runs are younger than the rest.
`StigmergyCommit-PT-11b` and `TokenRing-PT-050` were absent from `INPUTS/` --
they exceed the GitHub pages file size limit and `pnmcc-models-2026` does not
publish them -- so every one of their runs died on
`Cannot open file .../model.pnml`. The two archives were fetched from
`https://mcc.lip6.fr/2026/archives/INPUTS-2026.tar.gz`, which carries them at
full size, installed into the corpus and deployed, and their twelve tests rerun
against build `202609051741`. `TokenRing-PT-050` answers OneSafe, StableMarking
and all sixteen UpperBounds in under 15 s; `StigmergyCommit-PT-11b`, a 1.4 GB
PNML, answers OneSafe in 112 s and 8 of 16 bounds, and times out on the rest.

`ReachabilityCardinality`, `ReachabilityFireability`, the CTL and LTL
examinations and `StateSpace` were not run.

The `ReachabilityDeadlock` numbers below are those of this campaign. The whole
examination was rerun afterwards on a newer build; its logs live in
`/data/ythierry/MCC26archive/RD-2026-09-05/` and the rerun has a campaign
folder of its own.

## Verdicts against the 2026 oracle

| examination | oracle | known | answered | ok | wrong | missed | bonus | wall (h) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Liveness | 1953 | 1837 | 1813 | 1807 | 0 | 30 | 6 | 84.5 |
| OneSafe | 1953 | 1943 | 1942 | 1942 | 0 | 1 | 0 | 6.4 |
| QuasiLiveness | 1953 | 1783 | 1786 | 1763 | 0 | 20 | 23 | 128.9 |
| ReachabilityDeadlock | 1953 | 1893 | 1871 | 1868 | 0 | 25 | 3 | 44.6 |
| StableMarking | 1953 | 1886 | 1838 | 1837 | 0 | 49 | 1 | 65.5 |
| UpperBounds | 31248 | 30345 | 29982 | 29926 | 0 | 419 | 56 | 123.2 |

**No wrong verdict anywhere**: 39687 values compared, 0 disagreements with a
known oracle value. 89 answers are `bonus`, that is verdicts on formulas the
2026 consensus oracle leaves at `?` — candidate contributions to the oracle,
and the only answers no other tool checks, so the only place an error could
hide. They are the `status == bonus` rows of `verdicts.csv`.

Two different things are called a miss. 1781 values went unanswered — an ideal
tool would have produced every one of them. Of those, 1237 are `none`: nobody
answered, the oracle has `?` there too, and they are open problems rather than
a gap of ours. The other 544 are `missed`: the consensus knows the value and we
did not produce it. Only that second number measures us against the field.

| examination | oracle | answered | ok | missed | none | unanswered |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Liveness | 1953 | 1813 | 1807 | 30 | 110 | 140 |
| OneSafe | 1953 | 1942 | 1942 | 1 | 10 | 11 |
| QuasiLiveness | 1953 | 1786 | 1763 | 20 | 147 | 167 |
| ReachabilityDeadlock | 1953 | 1871 | 1868 | 25 | 57 | 82 |
| StableMarking | 1953 | 1838 | 1837 | 49 | 66 | 115 |
| UpperBounds | 31248 | 29982 | 29926 | 419 | 847 | 1266 |

`toolsupport.py` joins `verdicts.csv` with the contest's
`raw-result-analysis.csv` and says who did answer each of the 544. Every one is
backed by at least one tool and every `bonus` by none, so the oracle is fully
explained by the evidence. **307 of the 544 were answered by ITS-Tools itself
in the contest** — our run is weaker than the submission there, most plausibly
because the contest budget is larger than the 1800 s used here. The remaining
237 are values another tool has and ITS-Tools did not: Tapaal alone carries
most of them, and 43 come from `2025-gold`, the previous edition's medallist,
the only tool to answer them. Since this campaign the same information is in
the oracle itself: `pnmcc-models-2026` now writes
`TECHNIQUES ORACLE2026 ITSTOOLS TAPAAL 2025GOLD`, the tools whose own answer is
the consensus value.

## Why 544 values were not produced

| cause | L | OS | QL | RD | SM | UB | total |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| timeout (1800 s) | 28 | 0 | 17 | 17 | 15 | 338 | 415 |
| `petri64` not runnable | 0 | 0 | 0 | 0 | 20 | 20 | 40 |
| `OverlargeMarkingException` | 1 | 0 | 1 | 1 | 1 | 25 | 29 |
| eclipse fatal error | 0 | 0 | 0 | 7 | 5 | 3 | 15 |
| gave up before the timeout | 0 | 1 | 2 | 0 | 8 | 32 | 43 |
| `its-*` C++ abort | 1 | 0 | 0 | 0 | 0 | 0 | 1 |
| Java `OutOfMemoryError` | 0 | 0 | 0 | 0 | 0 | 1 | 1 |

Only the first and the last lines are about the tools being too slow or too
hungry; 85 of the 544 are harness or front-end accidents.

* **`petri64` not runnable.** Two windows, 15:32–15:56 and 19:09–19:24, where
  the Java side reports `Cannot run program .../bin/petri64: Exec failed`,
  first `No such file or directory` (1483 calls) and later `Permission denied`
  (146 calls). The binary in
  `itstools/itstools/plugins/fr.lip6.petrispot.binaries_*/bin/` was replaced
  while the campaign was running and came back mode `644`, and a run that
  cannot spawn `petri64` **loses PetriSpot silently** — the Java side catches
  the exception and carries on. A fresh install of the
  ITS-Tools product unzips **every** plugin binary as mode `644` — `petri64`,
  `its-reach-linux64`, `louvain-linux64` alike — so the runtime is the only
  thing that ever makes them executable, and under a hundred concurrent jobs
  sharing one tree that is a race. The redeploy of build `202609051741`
  `chmod +x`es the plugin `bin/` directories before the rsync, and a `p2.inf`
  in each binaries plugin of ITS-Tools now restores mode 755 at install time,
  which is the only point of the packaging chain that sees a real filesystem.
* **`OverlargeMarkingException`, libDDD overflow.** `GPPP-PT-C0010N1000000000`,
  `SharedMemory-COL-050000`, `SatelliteMemory-PT-X65535Y2048`: markings that do
  not fit the front-end's or libDDD's value type. An ITS-Tools limit.
* **Eclipse fatal error** on `HealthRecord-PT-15/16`, `HexagonalGrid-PT-126`,
  `MAPK-PT-10240` and 58 more logs. The referenced
  `itstools/configuration/*.log` files no longer exist on the cluster, so the
  cause is not recorded; re-running one of these instances with the
  configuration directory preserved is the way to find out.
* **Gave up before the timeout.** The run ends with `Detected timeout of ITS
  tools` or an `its-reach` error well before 1800 s and reports nothing
  (`BlocksWorld-PT-09..12` in StableMarking at 830–1300 s,
  `RERS17pb113-PT-7/9` in UpperBounds at ~1430 s). The per-solver budgets add
  up to less than the total budget: time is left on the table.

## What PetriSpot did

| examination | logs walking | walker calls | properties asked | solved | rate | walker (h) | killed | invariant timeouts |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Liveness | 647 | 1313 | 554330 | 181342 | 0.33 | 5.0 | 2 | 2 |
| OneSafe | 120 | 142 | 1397600 | 6149 | 0.004 | 1.2 | 0 | 1 |
| QuasiLiveness | 1389 | 2129 | 10814663 | 2175398 | 0.20 | 32.2 | 19 | 1 |
| ReachabilityDeadlock | 1434 | 1845 | 1845 | 929 | 0.50 | 0.4 | 0 | 0 |
| StableMarking | 1472 | 2430 | 3836373 | 2265340 | 0.59 | 27.0 | 19 | 0 |
| UpperBounds | 1919 | 10087 | 45851 | 19317 | 0.42 | 24.3 | 261 | 134 |

Read from the `PetriSpot walker: <solved>/<asked> properties solved in <t> ms`
lines. A "property" is atomic — one transition for QuasiLiveness, one place
pair for StableMarking, one place for OneSafe — so the counts are not
comparable across examinations, only the rates are.

90 h of walking out of 453 h of campaign: **one fifth of the whole compute
budget goes through the walk engine.** The walker is killed by its 50 s budget
in 301 calls, 261 of them in UpperBounds, where it is called five times per log
on average.

`PetriSpot timed out after 60 s (PSEMIFLOWS)` in 127 calls and `(PFLOWS)` in
11, concentrated on one family: `PolyORBLF-*` and `PolyORBNT-*` account for
101 of the 138, then `SafeBus-PT-1x/2x`, `FamilyReunion-*`,
`SharedMemory-COL-*`, `DatabaseWithMutex-PT-*`. That is the semi-flow solver's
own hard set on this benchmark.

## Where the walk engine could earn verdicts

The 415 timeouts are the honest target list. UpperBounds carries 338 of them
over 94 instances, and the families repeat: `MultiCrashLeafsetExtension-PT-*`,
`FunctionPointer-PT-*`, `DoubleExponent-PT-*`, `FamilyReunion-*`,
`DatabaseWithMutex-*`, `HouseConstruction-PT-*`, `Murphy-PT-*`,
`RERS17pb11x-PT-*`. In StableMarking and QuasiLiveness the same
`MultiCrashLeafsetExtension` and `RERS2020` families dominate.

On `RERS17pb113-PT-7` the shape of the failure is legible in the log: the walk
raises `Max Seen` from `[4,1,2,...]` to `[4,2,2,...]` and never touches
`Max Struct`, stuck at 7 for every expression, so no bound closes and 11 of 16
values go unreported after 92 s of walking. On `BlocksWorld-PT-10` the walker
spends 280 s for 0 of 53 properties.


## Log volume

419 MB of logs for 11706 runs, and 60 % of it is one printing habit. Counting
bytes by line shape:

| MB | lines | shape |
| ---: | ---: | --- |
| 110.8 | 915 029 | `Model ,\|S\| ,Time ,Mem(kb) ,fin. SDD ,...` |
| 69.6 | 1 468 890 | `Reachability property qltransition_N is true.` |
| ~71 | ~800 000 | the statistics rows themselves |
| 5.3 | 121 409 | `Invariant property smplace_N does not hold.` |
| 4.8 | 73 626 | `SDD proceeding with computation, N properties remain` |
| 3.7 | 14 331 | `Running PetriSpot : '<absolute path>' '--loadKERS=...'` |
| 3.5 | 135 589 | `Problem TDEADN is UNSAT` |

The first line is the CSV **header** of the statistics table, reprinted before
every single row: `Statistic::print_table` (libDDD `ddd/statistic.cpp`) is
header + line + trailer, and libITS calls it once per property
(`bin/main.cpp:642`, and again at `:588`). On `ErlangenMainframeV1-PT-bP11C06`
QuasiLiveness that is 69 273 headers in one 17.8 MB log. Printing the header
once per run and `print_line` per property would drop 110 MB of the 419 with no
loss of information — though `analysis/logs2csv.pl` matches the header to find
the line after it, so it would have to be adjusted.

The `is true.` line comes from `EarlyBreakObserver::update`
(libITS `bin/EarlyBreakObserver.hh:100`), one per property discharged, which is
legitimate but pays 48 bytes for a name we already have in the `FORMULA` line.

PetriSpot's own output is 11 MB of the 419, 2.6 %, about five lines per
invocation, and nothing in it repeats. The one avoidable line on our side is
the Java runner's `Running PetriSpot : <full command line>` (3.7 MB), which is
immediately followed by `[INFO] Running PetriSpot with arguments : [...]`
saying the same thing.
