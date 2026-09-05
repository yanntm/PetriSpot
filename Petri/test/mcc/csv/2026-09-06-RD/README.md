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

## What the difference is not

The two builds are three commits apart: a whitespace commit, the `p2.inf`
packaging fix and the `ITSRunner` echo gate. Neither of the latter two touches
a search: one changes file permissions, the other guards a `println`. The six
commits that separate `202609051237` from `202609051349` are the OneSafe SMT
budget, the OneSafe place naming, a parse error report and log line trimming.
Nothing in that set is the deadlock search.

So the exchange is not a deliberate behavioural change, and the most economical
reading is **run to run variance of a randomised search sitting at the edge of
its budget** — which the concentration of the gains in one family and the
losses in another argues against, so the reading is not settled. A repeat run
of the eight regressed instances on this same build is in `RDvar/` on the
cluster: if they answer this time, it is variance.

Either way there is a methodological consequence for the overhaul: **a single
campaign cannot measure a change worth less than about 1 % of verdicts**, and
this exchange is 35 instances out of 1953. A/B work needs repeats.

## The verbosity change

Confirmed active and irrelevant here: the statistics table appears 23 times
across the new campaign against 50 before, over 101 and 120 `its-*`
invocations. ReachabilityDeadlock never carried the noise the change targets —
zero `SDD proceeding with computation` lines and zero per property
announcements in either campaign — so the directory is 20 MB either way. The
examination that would show it is QuasiLiveness, whose 271 MB is 60 % those two
shapes, and it has not been rerun.
