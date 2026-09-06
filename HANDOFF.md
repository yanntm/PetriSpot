# Handoff — MCC 2026 campaign, the portfolio design, the first fixes, the total examinations

State of the work as of 2026-09-06, 11:30. Read `PORTFOLIO.md` first if you are
picking up the design, `TOTAL_QUERIES.md` for the total examinations; read
this for where everything lives and what is in flight.

## 2026-09-06, 21:30: the QLA rerun read, and what it taught

**QLA rerun collected** into `/data/ythierry/MCC26run/2026-09-06c/QLA`
(1 681 logs, 54 still running at 20:40; tables in
`Petri/test/mcc/csv/2026-09-06c/`, pages rebuilt with the set
`ITS-Tools latest` on `2026-09-06c/*` and the morning as
`ITS-Tools 2026-09-06`). On finished runs: 113 models better, 93 worse,
answered atoms 7.96 M -> 7.14 M. Gains where the quests shine
(ResIsolation 153 -> 49 754 per instance, ErlangenMainframe V0 to V2,
RERS17pb113); losses of two kinds, both reproduced locally and fixed on
master after product `202609061624`:

1. **The quest sweep loses where rarity won.** RERS17pb114-PT-1: 43 251
   atoms in 9 s this morning (the `==` goals were unquestable, the threads
   ran on rarity), 508 in the rerun (locally, 10 s: quests 4 088, rarity
   30 121). Fix 78160aa: `auto` spreads sync and rare over the threads with
   a restart policy per strategy (RERS 24 852; ResIsolation 1 000 targets in
   13.5 s instead of 6 to 8; Erlangen full QLA 1 564 instead of 2 792). A
   hedge: the design that replaces it is WALK_PLAN.md 10.11 (time sharing of
   many workers over few cores, shares that follow results), to comment.
2. **Calls of 100 to 340 s that solved nothing** (FamilyReunion, DLCflexbar,
   DLCshifumi, MultiCrashLeafsetExtension, ServersAndClients). The
   components built before the sweep: DLCflexbar-PT-8b's 3 040 semiflows
   share hub places and the edge tables took 41 s and 32 GB (the cluster's
   16 GB); FamilyReunion's phase 1 gives 113 229 flows and phase 2 never
   ends, phase 1 having no deadline. Fix 6a5d288: a covering of the smallest
   semiflows, refusal over `MAX_WORK`, phase 1 stops at the deadline with no
   basis, the sweep gets the total minus the components' time. Locally:
   DLCflexbar 76 160 of 76 160 in 35 s (8 979 in the rerun), FamilyReunion
   492 694 of 508 489 (45 444).

**The rerun in progress (SMA, UBA, RD, QL, SM, OS) runs product
`202609061624`, which carries both problems.** Its results on those families
are not the walker's fault to read; the chain (push done, `Inv-Linux`,
ITS-Tools rebuild, reinstall, rsync) must run again for master `6a5d288`
before the next campaign. Local models for the cases: `/data/ythierry/rers/`
(RERS17pb114-PT-1), `/data/ythierry/dlc/` (DLCflexbar-PT-8b),
`/data/ythierry/family/` (props for the example net), `/data/ythierry/dbmutex/`
and `/data/ythierry/philo/` (deadlock tests of the LP).

**LP engine** (`Petri/src/lp/`, commits 22b1d3c to a2e8ee9): sparse
product-form basis, rows as a `MatrixCol`, branches as row overlays,
`DeadlockRefiner`; DatabaseWithMutex PT-02 and PT-04 proved deadlock-free,
Philosophers up to 100 give the deadlock's Parikh vector, 1 000 and beyond
wait for the dual warm start. `--lpDebug` checks the rebuilt inverse.

## 2026-09-06, 20:00: the in-house LP engine, first prototype

`Petri/src/lp/` (design in its `algorithm.md`, map in `README.md`; commit
22b1d3c): a bounded-variable primal simplex over the state equation, goals
in disjunctive normal form, `--lp --lpHints=FILE --lpTime=S` in
`cli/LpDriver.h`. Tested as raw PetriSpot: Airplane's invariant proved;
49 of 50 ResIsolation QLA targets get a Parikh vector in about a second
(`Petri/test/logs/lp-resiso-50b.log`); Erlangen's three open targets in 2 s;
the Stigmergy target's hint takes the walk from 1.6 to 10 s to 1.2 s; the
CAN gathering target gets a vector the walk does not realise (spurious:
traps are the next refiner). Not built yet, in the order the design gives:
the exact rational check of an infeasibility certificate (until then the
`FORMULA` lines of `--lp` are a floating-point verdict for tests only, no
other path consumes them), the trap fixpoint and cut, a warm start of the
base problem across the atoms of a total examination, the dual simplex for
cuts, the rounding repair of a Parikh vector.

## 2026-09-06, 18:51: the rerun is submitted

Product `202609061624` ships the `petri64` of PetriSpot `0ae6e4c` (checksum
`36f11995e952`, matched against `Inv-Linux` and inside the deploy tree), with
the afternoon's walker work: `==` goals quested (7652b2e), the sweep restarts
on hopeless markings and evicts condemned pool starts (5ade6af), the stage
choice bounded and tabled (cda89d9), seeds fixed (929b720). Smoke through
the harness on ResIsolation-PT-N10P4 QLA, 90 s: 50 702 QLIVE verdicts.
Rsynced to the cluster at 18:50 while about 120 LTLF jobs were still
waiting (the user accepted the risk to those runs). The submitter
`~/MCC26/MCC-drivers/submit-2026-09-06c.sh`, detached on the cluster head,
log `submit-2026-09-06c.log`, runs `QLA SMA UBA RD QL SM OS` in that order,
1800 s / 4 cores / `tall`, one `run_oar.sh` after another; it moved the
previous results aside as `*-2026-09-06a` (QLA, SMA, UBA, RD from this
morning; QL, SM, OS from 2026-09-05). Those seven, then RC, RF, LTLC, L and
UB, were copied whole (stdout and stderr, file counts checked) to
`/data/ythierry/MCC26archive/{2026-09-06a,2026-09-06-baseline,2026-09-05}/`
and removed from the cluster, whose tree is the 29 GB of INPUTS plus the
running LTLF and the rerun; the stdout mirrors the pages read stay under
`/data/ythierry/MCC26run/`. LTLF is to be archived the same way when it
ends. L is not rerun. Collect into `/data/ythierry/MCC26run/2026-09-06c/` with the
rsync recipe and point `ITS-Tools latest` at it in
`~/git/MCC-analysis/campaign/example.json` (the 9-run warmup sets of the
09-05 totals were deleted; `build.py` no longer lists a total examination
twice). Yardsticks of the walker after today, 4 threads, local:
ResIsolation 1 000 targets claimed in 6 to 8.5 s; Erlangen full
QuasiLiveness 2 792 in 30 s (campaign 321); Stigmergy 1 496 (campaign 797).
Open design: WALK_PLAN.md 10.10, subgoals under a budget.

## 2026-09-06, 17:30: the product carries master, one walker fix on top

The CI chain ran once more: ITS-Tools product `202609061502` ships the
`petri64` of PetriSpot `42ad9cc` (checksum matched against `Inv-Linux`), with
both Java fixes; it is installed in `/data/ythierry/MCC26deploy/MCC-drivers/
itstools/` (install log `/data/ythierry/MCC26deploy/install-2026-09-06c.log`),
not yet rsynced to the cluster. LTLF was still draining (845 waiting at
17:00, about 450 jobs an hour); CTLC and CTLF were deleted, their directories
are empty.

**Found while smoking the product on ResIsolation-PT-N10P4 QLA** (harness,
90 s): 148 QLIVE verdicts, the campaign's level, while the same binary claims
50 000 targets standalone in 30 s. The Java hands the walker the atoms as
`(and (== p 1) ...)` over the pre-places (fireability on a safe net, 95 602
distinct pre-sets for 147 855 transitions), and `QuestSweep::distance` ranked
only `>=` atoms: every goal was UNREACHED, the threads walked on the rarity
filler. Fixed in PetriSpot `7652b2e` (`ComponentStrategy::questNeed` shared by
the quests and the ranking): 50 319 verdicts in the same 90 s. The captured
inputs are in `/data/ythierry/resiso/javahand/` (the `.pnet` and `.sexpr` the
Java wrote); logs `Petri/test/logs/yardstick-javahand-*.log`,
`/data/ythierry/MCC26deploy/smoke-resiso-QLA-*.log`. **The chain must run
again for `7652b2e`** (PetriSpot push, `Inv-Linux`, an empty commit on
ITS-Tools master, reinstall) before the rsync and the QLA rerun.

Also measured: the 1 000-target yardstick at `42ad9cc` claims 668 in 20 s
where the `a1b1732` binary claims 1 000 in 4.3 s (`Petri/test/logs/
yardstick-1000-*.log`); the regression is in `9d8cf28`..`42ad9cc`, untouched.

## State at the end of 2026-09-06 (about 16:45)

**Code.** PetriSpot master carries, beyond what the morning shipped: quests
that gather (k nearest tokens over the components holding the place), bounds
as a staircase in the quest sweep (`QuestSweep.h`, a bound quested at
known + 1), recursive staging (`ComponentStrategy.h` is a stack of stages,
see `Petri/src/walk/algorithm.md` and WALK_PLAN.md section 10), the pool never
offered a dead marking, hopeless targets left out of the sweep's ranking.
Controls pass (three ResIsolation joins 2.8 s, three Erlangen targets under
2 s, Stigmergy 1.7 s, Airplane); the 1 000-target sweep regressed from 1 000
in 7 s to about 800 in 30 s and is the first thing to recover (the flat
version is commit `3a379ec`'s `ComponentStrategy.h`); CANInsertWithFailure
PT-010's gathering is the open yardstick (`/data/ythierry/can/`, files
`ub-open.sexpr`, `gather.sexpr`; PT-005 is trivial). All of it is pushed.

**Cluster.** The CTL baseline jobs are deleted (3 906); about 1 000 LTLF jobs
remain queued. The rerun script `~/MCC26/MCC-drivers/submit-2026-09-06c.sh`
is staged, not run. The ITS-Tools product to deploy must ship a `petri64`
built from this master (the CI chain: PetriSpot push, `Inv-Linux` deploy,
an empty commit on ITS-Tools master, the product on lip6.github.io; check the
binary inside the zip for the `--runTime` help text before deploying). Then
`install_itstools.sh`, the rsync of `itstools/`, a warmup on a handful of
AirplaneLD jobs watched to the end, then the examinations one `run_oar.sh`
at a time. The user's rules: small batches first, never `xargs` a thousand
ids at `oardel`, watch before committing hundreds of CPU hours.

## Rerun in preparation (2026-09-06, 16:00)

The walker gained components and the quest sweep (PetriSpot up to `c514d28`),
ITS-Tools its two fixes (walker verdict polarity `3aeb7f95`, skeleton enabling
polarity `d26cddfa`). The chain to the cluster: PetriSpot CI deploys `petri64`
to `Inv-Linux` (done for `be8b012` at 13:43 UTC; `c514d28`, the boxed
semiflows, in progress), then ITS-Tools CI must rebuild so the product ships
that binary (an empty commit on ITS-Tools master triggers it; `cf1e02af` was
pushed for `be8b012` and is superseded; push another once `c514d28` is on
`Inv-Linux`), then `install_itstools.sh` in `/data/ythierry/MCC26deploy/
MCC-drivers` and the rsync of `itstools/` to the cluster (recipes below), a
warmup on AirplaneLD, then `~/MCC26/MCC-drivers/submit-2026-09-06c.sh` on the
cluster head (staged, not run: it moves `QLA SMA UBA RC RF RD` aside as
`*-2026-09-06a`, all mirrored locally, and resubmits them; edit `EXAMS`). The
queue held 5 500 baseline jobs (LTLF, CTLC, CTLF) at 15:45: the user decides
whether the rerun waits behind them. On the cluster the Java hands the QLA
atoms to the walker with a 30 s sweep, so the auto sweep choice (components,
quests) applies without a Java change. Collect the rerun into
`/data/ythierry/MCC26run/2026-09-06c/` and point `ITS-Tools latest` at it in
`MCC-analysis/campaign/example.json`.

## Next session: the walk engine on very large nets

Read `WALK_PLAN.md` sections 9 and 10 first. Section 9 measured why the walk
crawls on hub-dense nets (65 000 arc visits per step, the target index) and
what landed: per-walker up/down target lists, `--partition` (off), no
re-check after a reset to the initial marking, `distinct transitions fired`
in the thread report. Section 10 is the plan the user is to comment on
before code: SAT-style principles (restarts, counters, one adaptive
strategy, learning from blockages), the memory, the structure the invariant
engine gives (components from unit semiflows, interaction hypergraph,
projections, siphons), and eight ordered steps with yardsticks. Steps 1 and 2 landed
on 2026-09-06 (`Knowledge.h`, `NoveltyTracker.h`, `RestartPolicy.h`,
`RareStrategy.h`, commit 76bcd0a): restarts happen, the pool is fed, and the
yardstick moved from about 90 to 153 distinct transitions per thread, where
every 147 000-transition ResIsolation instance stops whatever the choice. The
family has no dead transition (N08P1 fully live, the field says quasi-live
wherever it decided), the net has 14 P-flows sharing p0..p12 and 13
token-producing transitions: the fourteen-input transitions are barriers.
**Step 4 landed the same day** (`Components.h`, `ComponentStrategy.h`,
`QuestSweep.h`, `--strategy=sync`, `--sweepChoice=sync`, commits 70e1a82 to
3a379ec): processes from the P-semiflows computed in-process, quests with
freezing, stages through the barriers whose outcome brings the goal closest,
a tabu on those that did not help, and a sweep that picks the nearest open
targets one per thread. Yardsticks at 30 s on 4 threads: the 1 000-target
file claimed in full in 7 s; the full QuasiLivenessAll of
ResIsolation-PT-N10P4 50 917 of 147 855 (campaign: 153 in 1800 s) and of
ErlangenMainframeV2-PT-bP10C08 1 570 of 79 094 (campaign: 321). The
`--debugSteps` trace prints the stage sequence and the stuck states. Next:
`sync` as an arm of the guided portfolio, the default sweep choice when
components exist, the other families of WALK_PLAN.md section 9, then the
self-configuration of section 10. The yardstick run:

```
cd /data/ythierry/resiso   # ResIsolation-PT-N10P4 extracted, qla-1000.sexpr from Petri/test/probes/qla_props.py
~/git/PetriSpot/build/petri64 -i ResIsolation-PT-N10P4/model.pnml --props=qla-1000.sexpr \
    --threads=4 --totalTime=20 --sweepTime=20 -t 20      # today: about 90 distinct transitions, 125 claims
```

## Session of 2026-09-06, afternoon: the campaigns read

* **Cluster.** RD, QLA, SMA, UBA complete (every oracle has exactly one run;
  the AirplaneLD warmup duplicates are set aside in
  `/data/ythierry/MCC26run/2026-09-06/warmup/`). Of the baseline submission,
  RC is complete and collected; RF was at 1902/1953 and LTLC, LTLF, CTLC, CTLF
  had not started at 11:00 (about 7900 jobs queued). Fetch them with the
  rsync recipe below, collect with `mcclogs2csv.py` into
  `Petri/test/mcc/csv/2026-09-06-baseline/`, rerun `report.py` there.
  Logs are in `/data/ythierry/MCC26run/2026-09-06/{RD,QLA,SMA,UBA,RC}` with
  their `.stderr`; the Eclipse `configuration/*.log` in `eclipse-logs/`.
* **Tables and reports** committed under `Petri/test/mcc/csv/2026-09-06/`
  and `csv/2026-09-06-baseline/`, each with a `README.md` (the reading) and a
  generated `REPORT.md` (every table). `Petri/test/mcc/README.md` documents
  the scripts: `report.py`, `toolboard.py`, `totalcheck.py`, `ubacheck.py`,
  `eclipselogs.py`.
* **A wrong verdict, reproducible.** RC `NeoElection-PT-3` formula 00 and
  `PT-7` formula 05: FALSE against a four tool TRUE. PetriSpot's
  `WalkDriver.h` `makeTargets` prints `FORMULA propN FALSE TECHNIQUES
  TOPOLOGICAL TRIVIAL` when a goal folds to constant false (a sound
  `Simplify.h` rule: the predicate over an emptied syphon), and ITS-Tools
  `PetriSpotWalker.readVerdicts` takes every FORMULA line as a witness
  (`publishWalkerVerdict`: "EF holds, AG does not") without reading the value.
  `INTEROP.md` line 409 promises PetriSpot never emits a non-witness verdict.
  Local reproduction: `Petri/test/logs/neo3-RC-local.log` (its-tools of the
  deploy on `/data/ythierry/neotest/NeoElection-PT-3`, 40 s). **Fixed in
  ITS-Tools** (commit "PetriSpotWalker: a FORMULA value carries the
  polarity"): `Verdicts.found` is WITNESS (1) or ABSENT (-1), every consumer
  reads the polarity; the user's choice, PetriSpot's FALSE is trusted as a
  proof. `INTEROP.md` section 5 documents it. Verified on a local product
  built from the parent pom (`mvn -o install -DskipTests` in
  `fr.lip6.move.gal.parent`, 1:22 min, product tarball under
  `ITS-commandline/fr.lip6.move.gal.itscl.product/target/products/`,
  extracted to `/data/ythierry/itstools-local/`): both instances answer TRUE
  (`Petri/test/logs/neo{3,7}-RC-fixed.log`). Not pushed. The ITS-Tools tree
  also carries the user's own uncommitted edits (`PetriSpotRunner.java`
  DEBUG=2 and `--useQPlusBasis`, MANIFEST.MF churn): left alone.
  **No oracle is published until the CI product carries the fix and RC and
  the totals are rerun.**
* **Total oracles, staged only.** `pnmcc-models-2026/merge_total_oracles.py`
  (uncommitted in that repo) merged the run vectors into the skeletons:
  `/data/ythierry/MCC26run/2026-09-06/oracles-merged/`, 0 conflicts, QLA
  51 %, SMA 89 %, UBA 54 % filled. Cross-checks against the consensus:
  `totalcheck.csv` 1484 QLA + 1603 SMA confirmed, `ubacheck.csv` 22 933
  bounds confirmed, 0 contradictions either way.
* **Web report.** Prototype in `~/git/MCC-analysis/campaign/` (pushed):
  `build.py example.json` builds one page per examination into
  `/data/ythierry/MCC26run/pages/` from result sets (our log directories,
  contest tools), `serve.py` serves them and the logs on 127.0.0.1:8080 for
  an SSH tunnel. The consensus, formula names and backing tools come from the
  oracle archive published by pnmcc-models-2026 (gh-pages, fetched to
  `/data/ythierry/MCC26run/oracle-2026/oracle`); the set `ITS-Tools latest`
  is a glob on `/data/ythierry/MCC26run/2026-09-06/*`, so an rsync plus
  `build.py example.json` shows the newest data. RF (1953) is in; LTLC was
  at 804 and LTLF, CTLC, CTLF not started at 12:20. The total examinations
  have their pages (`totals.py`: family progression, completion against
  time, atom by atom agreement between sets, implied global verdict against
  the consensus); every 2026 tool is a set, and a set that answered nothing
  on an examination is left off that page. Next: flag variants as sets,
  whatever the user asks after browsing.
* **RF: one wrong value shared with the contest.** `LastZero-COL-N20`
  ReachabilityFireability formula 13: ITS-Tools says FALSE in the contest and
  in our run, smpt, Tapaal and 2025-gold say TRUE. Deterministic on the
  ITS-Tools side, not the walker: the verdict is tagged `CPN_APPROX`, the
  skeleton over-approximation of the coloured net (4 places, 3 transitions)
  ruling the fireability unreachable before the unfolding (log
  `2026-09-06/RF/OAR.1336369`, line 106). **Fixed in ITS-Tools `d26cddfa`**
  (pushed): `Simplifier.allEnablingsAreNegated` asked for negated enablings
  whatever the root operator; the skeleton preserves enablings upward only,
  so an EF refuted there needs positive enablings and an AG proved there
  negated ones. Verified on a local product: LastZero-COL-N20, Philosophers-
  COL-000010, SharedMemory-COL-000005 RF all 16/16 against the oracle
  (`Petri/test/logs/*-RF-fixed.log`). In the RF campaign the bug bit once
  (18 skeleton verdicts, 17 sound AG TRUE, 1 EF FALSE wrong); RC's 2315
  skeleton verdicts are cardinality atoms, exact on the skeleton.
* **Findings to act on** (details in the two READMEs): UBA has 313 wall runs
  that never reach a walk, stalled after the invariants; QLA's residue is
  seven families where one SMT call on 79 k atoms eats 687 s for nothing;
  the contest ITS-Tools beats us on 40 RC instances, a dozen of them closed
  there in under three minutes while we burn 1800 s; `eclipse_fatal` is
  always a Java OutOfMemoryError (`-Xmx16384m`).

## Session of 2026-09-06: total examinations and streaming

* **ITS-Tools** gained `QuasiLivenessAll`, `StableMarkingAll`,
  `UpperBoundsAll` (`solver/total/`, P/T only), after a refactor of the global
  solver (`GlobalAtoms`, `ExhaustiveEngines`, `Aggregation`). Output is one
  line per atom, `QLIVE t12 TRUE <techniques>`, `BOUND p3 ? lo hi` while open.
* **Streaming.** PetriSpot walker verdicts are published into `DoneProperties`
  as the binary prints them (`PetriSpotWalker.Listener`), and PetriSpot prints
  a `BOUND` line each time a walk raises a bound. Before, everything a walk
  found was lost when the harness killed the run: Peterson-PT-5 bounds went
  from 17 to 396 of 834 closed in the same 120 s.
* **Oracles** of the total examinations are vectors, `T`/`F`/`?` per object
  wrapped at 80 columns (`pnmcc-models-2026/make_total_oracles.sh`, in
  `oracle.tar.gz`; 5043 files committed in MCC-drivers). `run_test.pl` reads
  them by index. `SupportedExamination.txt` in ITS-Tools-MCC admits the three.
* **Collectors** `Petri/test/mcc/totallogs2csv.py` (one row per run) and
  `log2oracle.py` (one run into its vector oracle).
* **Cluster.** Old `RD/` and the first total warmup archived in
  `/data/ythierry/MCC26archive/2026-09-06-precampaign/`, cleared on the
  cluster. Product `202609060003` deployed. AirplaneLD warmup for RD, QLA,
  SMA, UBA passed (RD 18/18 against the oracle, every atom closed on the
  totals). The four full campaigns, RD 1953 jobs then QLA, SMA, UBA 1681
  each, 1800 s / 4 cores / `tall`, were submitted from 02:20 on 2026-09-06 by
  `~/MCC26/MCC-drivers/submit-2026-09-06.sh` on the cluster head, sequential
  `run_oar.sh` calls, log `submit-2026-09-06.log` next to it. Results land in
  `RD/`, `QLA/`, `SMA/`, `UBA/`; read them with `mcclogs2csv.py` and
  `totallogs2csv.py`, warmup tables in
  `/data/ythierry/MCC26run/warmup-2026-09-06/csv/`. That submission ended at
  03:02; a second one, `submit-2026-09-06b.sh` (log `submit-2026-09-06b.log`),
  started at 03:32 with the baseline examinations RC, RF, LTLC, LTLF, CTLC,
  CTLF, 1953 jobs each, same settings, into `RC/` ... `CTLF/`. Read those with
  `mcclogs2csv.py` against `csv/2026-09-05/` (the campaign had no RC/RF/LTL/CTL
  then, so this is their first baseline).
* **RD, read early** (all 1953 in by 03:00, tables in
  `/data/ythierry/MCC26run/2026-09-06/csv-peek/`, baseline in
  `/data/ythierry/MCC26archive/2026-09-06-precampaign/csv/`): 0 wrong, missed
  14 -> 9, bonus 11 unchanged. Gained 8 deadlocks by the walk (the three Shield
  targets among them), lost 3 FALSE proofs (DatabaseWithMutex-PT-20, PGCD
  D02N100 PT and COL) to the wall: the escalated deadlock walk now spends 30 s
  per call, 4 to 6 calls per run, ahead of the engines that prove absence. The
  user's reading: walk beside the SMT proof, not before it, and check that the
  SMT deadlock timeouts escalate at all (`DeadlockTester`). Of the 9 missed,
  5 rest on a single tool (2 on 2025 gold, i.e. us last year).
* The `-timeout` flag of ITS-Tools is a per-engine budget, not a deadline: the
  harness kill ends every run that does not close its cohort.

### Picking up the 2026-09-06 campaigns

The submitter is detached on the cluster head (`setsid nohup`), the jobs are
OAR's: nothing depends on a local session. To pick up:

1. **Is the submission complete?** `tail ~/MCC26/MCC-drivers/submit-2026-09-06.log`
   ends with `SUBMISSION DONE` and four timestamps; `oarstat -u ythierry` says
   what is still queued or running. Expect RD 1953 logs, QLA, SMA, UBA 1681
   each. A job killed by walltime leaves a log without verdicts, a job never
   submitted leaves no log: diff the oracle list against the `Running test`
   lines as the `comm` recipe in `BENCH.md` does, and resubmit those alone.
2. **Fetch**, never with `--delete`:
   ```
   cd /data/ythierry/MCC26run && mkdir -p 2026-09-06
   for d in RD QLA SMA UBA ; do rsync -rz --exclude='*.stderr' cluster.lip6.fr:MCC26/MCC-drivers/$d 2026-09-06/ ; done
   ```
3. **Tables**, committed under `Petri/test/mcc/csv/2026-09-06/` with a
   `README.md` saying what the campaign showed:
   ```
   cd ~/git/PetriSpot/Petri/test/mcc
   python3 mcclogs2csv.py /data/ythierry/MCC26run/2026-09-06/RD -o csv/2026-09-06/
   python3 toolsupport.py ~/git/pnmcc-models-2026/website/raw-result-analysis.csv csv/2026-09-06/verdicts.csv -o csv/2026-09-06/support.csv
   python3 totallogs2csv.py /data/ythierry/MCC26run/2026-09-06/{QLA,SMA,UBA} -o csv/2026-09-06/ --oracles /data/ythierry/MCC26run/2026-09-06/oracles
   ```
   RD is read against `csv/2026-09-06-RD/` (the three Shield instances are the
   target). For the totals, `total-runs.csv` has per run: `completion`,
   `witnessed` against `proved`, the engine per verdict (`initial`, `walk`,
   `smt`, `dd`), the quartile times `t25..t100`, and `failure`. The questions
   `TOTAL_QUERIES.md` asks: where completion falls off along a family, whether
   the residue is the same atoms across instances (diff the vector oracles of
   two instances), knee or slope in the quartiles.

### Updating the oracles of pnmcc-models-2026 from our logs

Today every published vector is `?`. The pipeline to fill them:

1. `log2oracle.py -o DIR QLA/*.stdout` (or the collector's `--oracles`) writes
   one vector per run, the run's verdicts and `?` where it said nothing.
2. **Merge** into the published vectors: a `?` takes the run's value, an equal
   value stays, a disagreement (`T` against `F`, two different bounds) is a
   finding to report, never overwrite. This script does not exist yet; it is a
   few lines over two vectors of the same length, position by position, and
   belongs in `pnmcc-models-2026` next to `make_total_oracles.sh`.
3. **Store** the filled vectors in `pnmcc-models-2026`: the generator runs at
   CI time and produces `?` only, so a committed folder of filled files, copied
   over the generated ones before `oracle.tar.gz` is packed, is the shape (the
   same way `oracleSS.tar.gz` is layered in by `install_inputs.sh`). Then
   `MCC-drivers/oracle/` is refreshed from the archive and committed, as for
   the consensus oracles.

A vector filled from our own runs is a regression oracle, not a truth: the
self certifying side (QLIVE `T`, STABLE `F`, a bound's lower end) is
witnessed, the other side rests on our proofs alone. `TOTAL_QUERIES.md`
says so; keep that caveat with the files.

## What this session did

1. Collected the first MCC 2026 cluster campaign and built the tooling to read
   it (`Petri/test/mcc/`).
2. Repaired the harness: the two oversized models, the non-executable
   `petri64`, the log volume.
3. Wrote `PORTFOLIO.md` — a design for the solving loop with goals G1..G7 — and
   revised it twice against the code and against the user's answers.
4. Implemented the first two fixes: `--escalate` in PetriSpot, and a walk
   running beside the decision diagrams in ITS-Tools.

## Where things are

| what | where |
| --- | --- |
| campaign logs, current | `/data/ythierry/MCC26run/{L,OS,QL,RD,SM,UB}` |
| superseded RD logs, and the variance recheck | `/data/ythierry/MCC26archive/RD-2026-09-05/`, `.../RDvar/` |
| the harness, built locally then rsynced | `/data/ythierry/MCC26deploy/MCC-drivers/` |
| the cluster tree | `cluster.lip6.fr:~/MCC26/MCC-drivers/` |
| collected tables, committed | `Petri/test/mcc/csv/2026-09-05/`, `csv/2026-09-06-RD/` |
| Shield models extracted for local tests | `/data/ythierry/shieldtest/` |
| maven logs of this session | `/data/ythierry/itstools-mvn{5..10}.log` |

Repositories touched: `PetriSpot` (this one), `lip6/ITSTools`
(`~/git/ITStools`), `yanntm/pnmcc-models-2026`. `BENCH.md` documents the
campaign harness; `Petri/test/mcc/README.md` the collectors.

## Tools written this session

* `Petri/test/mcc/mcclogs2csv.py` — a directory of MCC logs into `runs.csv`
  (one row per log: verdict census, failure signature, PetriSpot and walker
  counters) and `verdicts.csv` (one row per formula). The oracle is *in* the
  log, so no oracle file is needed.
* `Petri/test/mcc/toolsupport.py` — joins `verdicts.csv` with the contest's
  `raw-result-analysis.csv` to say which tools produced each value we missed.
* `Petri/test/probes/perprop_payoff.py` — what the per property reduction phase
  of UpperBounds costs and returns.
* `Petri/test/probes/reduction_resistance.py` — first-pass reduction yield
  against whether the run answered everything.
* `Petri/test/probes/dd_tail.py` — how long a run spends past the Java side's
  last word, i.e. in the single threaded diagram engine.

## The numbers that drive the design

From `Petri/test/mcc/csv/2026-09-05/` (1953 instances, six examinations,
1800 s, 4 cores) unless said otherwise.

* **0 wrong verdicts** in 39 687 comparisons. 544 values the consensus has and
  we do not; 1237 nobody has; 89 we have and the consensus does not.
* The walker got **90 h of a 1812 core-hour budget, 5 %**.
* In the 601 runs that used the whole 1800 s, the tail past the Java side's
  last word is **48 % (Liveness), 41 % (QuasiLiveness), 31 % (StableMarking)**.
* The UpperBounds per property phase: 250 attempts, 68 productive,
  **47 319 s spent, 4 502 s of it productive**.
* Reduction yield separates the runs that answer everything from the others
  (38.9 % against 22.0 % of places removed on pass 1) but not cleanly enough to
  gate anything alone.
* 419 MB of logs, 60 % of it a statistics header reprinted per row and one
  verdict announcement per property.

## In flight

**ITS-Tools CI.** The last push (`4f392b34` plus the follow-up described below)
must produce a product before any rerun. `petri64` is fetched from the
PetriSpot `Inv-Linux` branch **at ITS-Tools build time**, so always push
PetriSpot first, wait for its CI to deploy, then push ITS-Tools. Otherwise a
new Java flag meets an old binary and CLI11 rejects it, killing every walk.
Check with:

```
git fetch origin Inv-Linux && git log -1 --oneline origin/Inv-Linux
curl -sIL https://lip6.github.io/ITSTools/fr.lip6.move.gal.itscl.product-linux.gtk.x86_64.zip | grep -i last-modified
```

**Cluster.** See the 2026-09-06 section above for what is queued.

## What to do next

1. **Rerun ReachabilityDeadlock** against `csv/2026-09-06-RD/` once the CI has
   published. The three consistently lost instances (`ShieldIIPt-PT-010B`,
   `ShieldIIPt-PT-020B`, `ShieldPPPt-PT-030B`) are the target; locally,
   `--escalate` takes `ShieldPPPt-PT-030B` from 0/3 to 3/3 and finds
   `ShieldIIPt-PT-010B` at 23.5 s on the third escalation.
2. **Watch UpperBounds and QuasiLiveness for a slowdown.** `--escalate` is now
   on for the reachability and bounds walks too, not only deadlock. It makes a
   walk spend its whole `--totalTime` instead of returning early, which is the
   point, but it also delays the outer loop's next reduction pass. If UB wall
   time explodes, that is why, and the flag to reconsider is the `--escalate`
   argument in `PetriSpotWalker.runReachability` / `runBounds`.
3. **Measure the log volume on QuasiLiveness.** The `ITSRunner` verbosity gate
   was confirmed active but RD never carried the noise; QL's 271 MB is where it
   should show.
4. **G1 of `PORTFOLIO.md`**: effort accounting in `DoneProperties`. Everything
   else in the plan is argued from measurements we do not yet collect from the
   tool itself.

## Fixes landed

**PetriSpot `c064fa8`** — `--escalate`. `WalkDriver.h` stopped the rounds when
one "solved nothing, improved no bound, and every walk ended on the step
budget". True of *time*, false of *steps*: the flag multiplies the step budget
and walks on while `--totalTime` lasts.

**ITS-Tools `847311e9`** — `runDeadlock` passes `--totalTime` and `--escalate`.
It was the only walker entry point passing neither, so the deadlock walk did
one round and handed the rest back: six calls of 1045, 643, 337, 1146, 409 and
371 ms against a 35 s budget, then 1796 s in `its-ctl`. `PetriSpotWalker` also
gained a `Cancel` handle, a per-call thread count and `runBeside`.

**ITS-Tools `4f392b34`** — `solver/ParallelWalk.java`: a three thread PetriSpot
walk beside every diagram attempt, publishing into the shared
`DoneProperties`. `ParallelWalk.ENABLED` turns it off. `verifyWithSDD` is split
into a wrapper that starts and stops it around `runExhaustiveEngines`, which
holds the former body and its several early exits.

It declines to start when `computeToCheck` cannot state every open property as
a target, rather than walk a partial set — so the deadlock examination has no
companion walk yet, its property not being an EF/AG shape. Giving it one means
routing that case to `runDeadlock`.

**ITS-Tools, same push** — `timeout -= (System.currentTimeMillis() - time) *
1000` in `runExhaustiveEngines` multiplied milliseconds by a thousand, so the
value went hugely negative and the `if (timeout < 0) return` right after it
meant **LTSmin was never reached**. Now `/ 1000`.

## Harness facts worth not rediscovering

* `StigmergyCommit-PT-11b` and `TokenRing-PT-050` are cut from the gh-pages
  corpus by a file size limit. Get them from
  `https://mcc.lip6.fr/2026/archives/INPUTS-2026.tar.gz` (1.1 GB, per-model
  `.tgz` inside). Both are now in `~/git/pnmcc-models-2026/website/INPUTS/` and
  on the cluster.
* **Never `rsync --delete` the harness root** once a campaign has run: the
  result directories exist only on the cluster. Rsync a subtree.
* A product built before ITS-Tools `202609052009` ships every plugin binary
  mode `644`; `chmod +x .../plugins/*/bin/*` before the rsync. From that build
  on, a `p2.inf` chmod touchpoint does it at install time.
* The oracle files of `pnmcc-models-2026` now name the tools that produced each
  consensus value (`TECHNIQUES ORACLE2026 ITSTOOLS TAPAAL 2025GOLD`).
* An `eclipse fatal error` hits this tree at random — 62 logs in the September 5
  campaign — and looks exactly like a tool failure in the CSV. Check the
  `failure` column before believing a one second miss.
* **Run to run variance is real.** Seven of sixteen repeats of eight marginal
  deadlock instances disagreed with the campaign. No A/B worth less than a few
  dozen instances can be read from a single campaign run.

## Useful commands

```
# the total examinations: tables, and one vector oracle per run
python3 ~/git/PetriSpot/Petri/test/mcc/totallogs2csv.py QLA SMA UBA -o <outdir> --oracles <outdir>/oracles
python3 ~/git/PetriSpot/Petri/test/mcc/log2oracle.py QLA/OAR.123.stdout

# collect a campaign into CSV, from the log directory
cd /data/ythierry/MCC26run
python3 ~/git/PetriSpot/Petri/test/mcc/mcclogs2csv.py L OS QL RD SM UB -o <outdir>
python3 ~/git/PetriSpot/Petri/test/mcc/toolsupport.py \
    ~/git/pnmcc-models-2026/website/raw-result-analysis.csv <outdir>/verdicts.csv -o <outdir>/support.csv

# fetch logs from the cluster (never with --delete)
rsync -rz --exclude='*.stderr' cluster.lip6.fr:MCC26/MCC-drivers/RD .

# deploy a fresh ITS-Tools, then push it to the cluster
cd /data/ythierry/MCC26deploy/MCC-drivers && rm -rf itstools && ./install_itstools.sh \
  && (cd itstools && ./install.sh)
rsync -rlptD --no-g --chmod=Dg+s --delete /data/ythierry/MCC26deploy/MCC-drivers/itstools/ \
  cluster.lip6.fr:MCC26/MCC-drivers/itstools/

# submit one examination (1953 jobs, from the foreground, let it finish)
ssh cluster.lip6.fr 'cd ~/MCC26/MCC-drivers && TIMEOUT=1800 WALLTIME=0:35:0 CORES=4 \
  HOSTS="tall%" ./run_oar.sh "oracle/*-RD.out"'

# a local deadlock walk, the fast loop used to test --escalate
cd /data/ythierry/shieldtest && echo '(deadlock prop0)' > dl.sexpr
~/git/PetriSpot/build/petri64 -i ShieldIIPt-PT-010B/model.pnml --props=dl.sexpr \
  --walkSteps=1250000 -t 30 --threads=4 -q --escalate --totalTime=30
```
