# Handoff — MCC 2026 campaign, the portfolio design, the first fixes, the total examinations

State of the work as of 2026-09-06, 11:30. Read `PORTFOLIO.md` first if you are
picking up the design, `TOTAL_QUERIES.md` for the total examinations; read
this for where everything lives and what is in flight.

## Next session: the walk engine on very large nets

Read `WALK_PLAN.md` sections 9 and 10 first. Section 9 measured why the walk
crawls on hub-dense nets (65 000 arc visits per step, the target index) and
what landed: per-walker up/down target lists, `--partition` (off), no
re-check after a reset to the initial marking, `distinct transitions fired`
in the thread report. Section 10 is the plan the user is to comment on
before code: SAT-style principles (restarts, counters, one adaptive
strategy, learning from blockages), the memory, the structure the invariant
engine gives (components from unit semiflows, interaction hypergraph,
projections, siphons), and eight ordered steps with yardsticks; step 0 and 1
are the natural start. The yardstick run:

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
