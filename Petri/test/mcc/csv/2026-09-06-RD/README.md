# ReachabilityDeadlock, rerun on build 202609052009

The whole examination rerun after the packaging fix (`p2.inf`, so `petri64`
ships executable) and the `ITSRunner` verbosity change, same settings as the
campaign of `../2026-09-05/`: 1953 instances, 1800 s, 4 cores, `tall` nodes.
The superseded logs are in `/data/ythierry/MCC26archive/RD-2026-09-05/`.

|  | old | new |
| --- | ---: | ---: |
| build | 202609051237 (99 %) | 202609052009 |
| known oracle values | 1893 | 1893 |
| answered | 1871 | 1890 |
| ok | 1868 | 1879 |
| wrong | 0 | 0 |
| missed | 25 | **14** |
| unanswered by anyone | 57 | 49 |
| bonus | 3 | **11** |
| runs where `petri64` could not be spawned | 1 | **0** |
| wall | 44.6 h | 36.9 h |
| walker | 0.45 h over 1845 calls | 1.12 h over 1755 calls |

Net: eleven fewer misses and eight more values than the 2026 consensus has,
in eight hours less wall time. Still no wrong verdict.

## The exchange behind the net figure

Thirty-five instances changed status, in both directions.

* **19 newly correct**, all from two families: the sixteen `ResIsolation-PT-*`
  and `Angiogenesis-PT-20/25/50`, plus `PGCD-COL-D02N100`. Every one of them
  was a timeout or a give-up before.
* **8 new bonus**, the whole `RERS17pb113-PT-2..9` run, values the consensus
  oracle leaves at `?`.
* **8 regressions**: `ShieldPPPt-PT-030B/040B/050B/100B`,
  `ShieldIIPt-PT-010B/020B`, `RERS2020-PT-pb104` and `ASLink-PT-07a`. The first
  seven answered in 3 to 161 s before and now burn the full 1800 s;
  `ASLink-PT-07a` dies in 1 s on an eclipse fatal error.

`ShieldPPPt-PT-040B` is the clean specimen. Both runs reduce identically
through the same rules; the old one then prints
`FORMULA ReachabilityDeadlock TRUE TECHNIQUES TOPOLOGICAL STRUCTURAL_REDUCTION
RANDOM_WALK` after a single walk call and four seconds. The new one makes six
walk calls totalling four seconds, finds nothing, and spends the remaining
1796 s inside `its-ctl` until the kill. The deadlock did not become harder to
reach: it stopped being found in the first sweep, and nothing in the run
reinvested in looking.

## What the difference is

Our own deadlock heuristic. `petri64` is fetched from the `Inv-Linux` branch
when the ITS-Tools CI builds, so the plugin carries whatever PetriSpot had
published at that moment:

| | plugin | petri64 published from |
| --- | --- | --- |
| old campaign | `202609051237` | before `e865bf5` |
| new campaign | `202609052009` | `45612fb`, deployed 14:50 UTC |

Between the two, three commits: `058dd89` designed the deadlock distance,
`e865bf5` added `walk/DeadlockStrategy.h` — a gradient toward a deadlock target
and saturating pools — and `45612fb` routed hints into the deadlock search and
made successor scoring degrade gracefully. The old binary had no
`DeadlockStrategy` at all; RD walks used the generic random and best-first
rounds.

So this rerun is an A/B of that work, and it reads: **+19 instances, +8 values
the consensus lacks, −8 instances**. The `RERS17pb113` family, which is hard,
went from nothing to eight answers.

## Why the eight came back empty

`ShieldPPPt-PT-040B` again. The new run calls the walker six times and each
call returns in **1045, 643, 337, 1146, 409 and 371 ms** having solved 0 of 1 —
against a budget of 35 s. It is not running out of time; it is stopping.
`WalkDriver.h` ends the rounds when one of them "solved nothing, improved no
bound, and every walk in it ended on the step budget", which on a single
property is reached in about a second. The old binary did not need the rule:
its first round found the deadlock in 426 ms.

Then ITS-Tools, having nothing from the walk, hands the model to `its-ctl`,
which runs 1796 s and is killed. Thirty-four of the walker's thirty-five
seconds, and all of the remaining half hour, went unspent on the search that
used to succeed here.

That is the same failure as the outer loop's no-progress exit, one level down:
a strategy that finds nothing quickly concludes rather than reinvests. It is
also the cheapest thing to try — the rule is one condition in `WalkDriver.h`.
A repeat of the eight instances on this same build is running in `RDvar/` on
the cluster to separate the heuristic from run-to-run variance.

## The verbosity change

Confirmed active and irrelevant here: the statistics table appears 23 times
across the new campaign against 50 before, over 101 and 120 `its-*`
invocations. ReachabilityDeadlock never carried the noise the change targets —
zero `SDD proceeding with computation` lines and zero per property
announcements in either campaign — so the directory is 20 MB either way. The
examination that would show it is QuasiLiveness, whose 271 MB is 60 % those two
shapes, and it has not been rerun.
