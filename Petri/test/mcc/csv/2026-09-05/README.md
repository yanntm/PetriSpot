# Campaign of 2026-09-05 — global properties and UpperBounds

ITS-Tools (`BK_TOOL=itstools`, no reducer prefix) on the whole MCC 2026 set,
1953 instances, six examinations, one `oarsub` per pair, 1800 s per test, 4
cores, all on the `tall` nodes of `cluster.lip6.fr`. Builds `202609051237`
(the first jobs) and `202609051349` (the rest). Logs live in
`/data/ythierry/MCC26run/`, rsynced from `~/MCC26/MCC-drivers/{L,OS,QL,RD,SM,UB}`.

`ReachabilityCardinality`, `ReachabilityFireability`, the CTL and LTL
examinations and `StateSpace` were not run.

## Verdicts against the 2026 oracle

| examination | oracle | known | answered | ok | wrong | missed | bonus | wall (h) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Liveness | 1953 | 1837 | 1813 | 1807 | 0 | 30 | 6 | 84.0 |
| OneSafe | 1953 | 1943 | 1940 | 1940 | 0 | 3 | 0 | 6.4 |
| QuasiLiveness | 1953 | 1783 | 1786 | 1763 | 0 | 20 | 23 | 127.9 |
| ReachabilityDeadlock | 1953 | 1893 | 1871 | 1868 | 0 | 25 | 3 | 43.6 |
| StableMarking | 1953 | 1886 | 1837 | 1836 | 0 | 50 | 1 | 65.4 |
| UpperBounds | 31248 | 30345 | 29958 | 29902 | 0 | 443 | 56 | 122.3 |

**No wrong verdict anywhere**: 39116 values compared, 0 disagreements with a
known oracle value. 89 answers are `bonus`, that is verdicts on formulas the
2026 consensus oracle leaves at `?` — candidate contributions to the oracle,
and the only answers no other tool checks, so the only place an error could
hide. They are the `status == bonus` rows of `verdicts.csv`.

## Why 571 values were not produced

| cause | L | OS | QL | RD | SM | UB | total |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| timeout (1800 s) | 28 | 0 | 16 | 16 | 15 | 330 | 405 |
| `petri64` not runnable | 0 | 0 | 0 | 0 | 20 | 20 | 40 |
| model absent from `INPUTS/` | 0 | 2 | 1 | 1 | 2 | 32 | 38 |
| `OverlargeMarkingException` | 1 | 0 | 1 | 1 | 1 | 25 | 29 |
| eclipse fatal error | 0 | 0 | 0 | 7 | 5 | 3 | 15 |
| gave up before the timeout | 0 | 1 | 2 | 0 | 7 | 32 | 42 |
| `its-*` C++ abort | 1 | 0 | 0 | 0 | 0 | 0 | 1 |
| Java `OutOfMemoryError` | 0 | 0 | 0 | 0 | 0 | 1 | 1 |

Only the first and the last lines are about the tools being too slow or too
hungry; 123 of the 571 are harness or front-end accidents.

* **`petri64` not runnable.** Two windows, 15:32–15:56 and 19:09–19:24, where
  the Java side reports `Cannot run program .../bin/petri64: Exec failed`,
  first `No such file or directory` (1483 calls) and later `Permission denied`
  (146 calls). The binary in
  `itstools/itstools/plugins/fr.lip6.petrispot.binaries_*/bin/` was replaced
  while the campaign was running and came back mode `644`. It is still `644`
  on the cluster, so **every future run there loses PetriSpot silently** — the
  Java side catches the exception and carries on. `chmod +x` before the next
  campaign, and check what `install.sh` does with the mode.
* **Model absent from `INPUTS/`.** `StigmergyCommit-PT-11b` and
  `TokenRing-PT-050`, the two archives missing from `pnmcc-models-2026`
  (BENCH.md). Known, unfixable from the models repository.
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

90 h of walking out of 449 h of campaign: **one fifth of the whole compute
budget goes through the walk engine.** The walker is killed by its 50 s budget
in 301 calls, 261 of them in UpperBounds, where it is called five times per log
on average.

`PetriSpot timed out after 60 s (PSEMIFLOWS)` in 127 calls and `(PFLOWS)` in
11, concentrated on one family: `PolyORBLF-*` and `PolyORBNT-*` account for
101 of the 138, then `SafeBus-PT-1x/2x`, `FamilyReunion-*`,
`SharedMemory-COL-*`, `DatabaseWithMutex-PT-*`. That is the semi-flow solver's
own hard set on this benchmark.

## Where the walk engine could earn verdicts

The 405 timeouts are the honest target list. UpperBounds carries 330 of them
over 93 instances, and the families repeat: `MultiCrashLeafsetExtension-PT-*`,
`FunctionPointer-PT-*`, `DoubleExponent-PT-*`, `FamilyReunion-*`,
`DatabaseWithMutex-*`, `HouseConstruction-PT-*`, `Murphy-PT-*`,
`RERS17pb11x-PT-*`. In StableMarking and QuasiLiveness the same
`MultiCrashLeafsetExtension` and `RERS2020` families dominate.

On `RERS17pb113-PT-7` the shape of the failure is legible in the log: the walk
raises `Max Seen` from `[4,1,2,...]` to `[4,2,2,...]` and never touches
`Max Struct`, stuck at 7 for every expression, so no bound closes and 11 of 16
values go unreported after 92 s of walking. On `BlocksWorld-PT-10` the walker
spends 280 s for 0 of 53 properties.
