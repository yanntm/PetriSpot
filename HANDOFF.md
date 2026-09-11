# Handoff — current state, next actions

Read this first, then the design file of the thread you pick up. This file is
rewritten, never appended to: what is done leaves it (result in the README or
the design file, history in git and `docs/HISTORY.md`). State as of
2026-09-11, evening.

## Orientation

| thread | read | code |
| --- | --- | --- |
| native structural reductions | `PS_REDUCTIONS.md`, `Petri/src/reduction/README.md`, `Petri/src/reduction/algorithm.md` | `Petri/src/reduction/`, `Petri/test/reduction/` |
| the walk engine (reachability) | `WALK_PLAN.md` sections 9, 10 | `Petri/src/walk/` |
| the solving loop, walker budgets | `PORTFOLIO.md`, `INTEROP.md` | `cli/WalkDriver.h`, ITS-Tools `PetriSpotWalker` |
| the CTL checker | `CTL_PLAN.md` sections 9, 11 | `Petri/src/ctl/`, `cli/CtlDriver.h` |
| the state equation, LP hints | `Petri/src/lp/algorithm.md` | `Petri/src/lp/` |
| the campaigns, the harness | `BENCH.md`, `Petri/test/mcc/README.md`, `TOTAL_QUERIES.md` | `Petri/test/mcc/` |
| libHSC as a competitor and companion | `HSC_PLAN.md`, `HSC_EXPERIMENTS.md`, `libHSC_in_MCC.md`; state in libHSC `handoff_mcc.md` | `~/git/MCC-drivers/hsc/`, ITS-Tools `hsc/`, `interop/` |
| operating guides: cluster, CI chain, local builds | `docs/CLUSTER.md`, `docs/CI.md`, `docs/BUILD.md` | |

Repositories: this one; `~/git/ITStools` (origin `lip6/ITStools`, the Java
side); `~/git/MCC-drivers` (the harness, origin `yanntm/MCC-drivers`, the
`itstools/` driver is a clone of `yanntm/ITS-Tools-MCC`); `~/git/MCC-analysis`
(the result pages); `~/git/pnmcc-models-2026` (models and oracles). The
ITS-Tools product bundles the `petri64` of the PetriSpot `Inv-Linux` branch
*at build time*: push PetriSpot first, wait for its CI, then push ITS-Tools.

## Where things are

| what | where |
| --- | --- |
| campaign logs, collected | `/data/ythierry/MCC26logs/<tool>/<build>/<EXAM>/`, `csv/` beside them; `local/<name>/` a local experiment's files |
| the deploy tree of the harness | `/data/ythierry/MCC26deploy/MCC-drivers/` (product `202609080208`) |
| the cluster tree | `cluster.lip6.fr:~/MCC26/MCC-drivers/`, results only there |
| result pages | `/data/ythierry/MCC26logs/web/campaign`, built by `~/git/MCC-analysis`; `web/order-sweep` the libHSC sweep's |
| native image material | `~/git/ITStools/ITS-commandline/native/`, `cluster.lip6.fr:MCC26/flat-test/`, the x86-64-v2 image `/data/ythierry/MCC26deploy/native-v2/` |
| GraalVM for a local image build | `/data/ythierry/graal/graalvm-jdk-25.0.4+7.1` |
| development models | `bench/models/<model>/` (git-ignored) |

## In flight

**Campaign itstools/202609101624: ITS-Tools on RC and RF, PT nets, small%,
complete.** Submitted 2026-09-11 11:30 from the cluster head by
`Petri/test/mcc/submit-2026-09-11.sh` (its log ends with `SUBMISSION DONE`):
1953 jobs an examination, 1800 s, 6 cores, `TAG=itstools`, the CI product
202609101624 through the locally built x86-64-v2 native image. Every job
has run (`RC.itstools` and `RF.itstools` hold 3906 files each, the queue is
empty). Not collected yet:

```
bash Petri/test/mcc/collect.sh itstools/202609101624 RC.itstools RF.itstools --pages
```

then the README of `MCC26logs/itstools/202609101624/`, its line in the
`MCC26logs/README.md` table, and a set for it in `campaign/example.json`
before `--pages` picks it up. Its `petri64` predates the preparation
pipeline (`reduction/Pipeline.h`): it is the baseline the next PetriSpot
campaign measures against, not a measure of the pipeline.

**Campaign hsc/20260910: libHSC alone on CTLC and CTLF, PT nets, small%.**
Submitted 2026-09-10 18:28 from the cluster head by
`Petri/test/mcc/submit-2026-09-10.sh` (detached, its log
`~/MCC26/MCC-drivers/submit-2026-09-10.log` ends with `SUBMISSION DONE`):
CTLC then CTLF over `oracle/*-PT-*-<EXAM>.out`, 1681 jobs each, 600 s,
6 cores, `HOSTS=small%`, `TAG=hsc`, results in `CTLC.hsc` and `CTLF.hsc`.
The tool is `hsc-pn` from libHSC `142fd4b`, the first CTL run since the
three defects the order sweep exposed were fixed (a stopped closure decides
nothing, a partial reachable set answers only what stands, a deadline is
never an error; libHSC `c20f91f`..`142fd4b`). Six cores because OAR caps a
job's memory at the node's RAM per core times the cores asked, and 6 of
small's 24 hyper-threads are 16 GB (`docs/CLUSTER.md` section 1). The
warmup (AirplaneLD-PT-0010, the eight examinations the tool declares) is in
`/data/ythierry/MCC26logs/hsc/20260910/_warmup/`: every one answered on
`small10`, no regression line, no failure. The cluster also holds
`CTLC.hsc600` and `CTLF.hsc600`, one Airplane warmup log each from the hsc600 campaign.

Watch, collect, read:

```
bash Petri/test/mcc/cluster_status.sh CTLC.hsc CTLF.hsc
bash Petri/test/mcc/collect.sh hsc/20260910 CTLC.hsc CTLF.hsc --pages
```

then the README of `MCC26logs/hsc/20260910/` and its line in the
`MCC26logs/README.md` table, and free the cluster (`docs/CLUSTER.md`
sections 4, 5). What to read: the wrong verdicts first (the sweep had 559 on
CTLC, 97 % off a partial set; the fix must bring that to 0), then how many
formulas the checker answers against the hsc600 and the ITS-Tools CTL
campaigns on the pages. A wrong verdict now is a new bug, not the old one.

**Campaign 202609080313, the first CTL examinations at scale (ITS-Tools,
tall%)**, is collected: `MCC26logs/itstools/202609080313/` (CTLC, CTLF, L,
5862 logs). Its reading is still to do, against the field (`report.py
--raw`, the pages against `ITS-Tools 2026` and `Tapaal 2026`): how many
formulas the explicit checker answers, how many only it had before the
diagrams, any wrong verdict. Its companion `hsc-pn` was the pre-fix binary,
so its `-hsc` verdicts on CTL carry the partial-set bug: a wrong one there
is explained before it is investigated.

**CI.** Both fixes are published and deployed in `202609080313`: `7f4113e0`
puts `--no-V` back in the LTSmin runners, `fcdae5b8` registers the coloured
`Sort[]`.

**The native image is now the launcher.** `runeclipse.sh` execs
`its-tools-native` whenever the file is present, so the deploy decides it;
`docs/CLUSTER.md` section 1 says how to tell and how to choose. The CI image
is AVX2, hence `tall%` only; a `NATIVE_MARCH=x86-64-v2` build runs on `small%`
and `big%` (`docs/CLUSTER.md` section 1, the node table). Its closed world was missing
`fr.lip6.move.gal.InstanceDecl[]` and `fr.lip6.move.gal.Synchronization[]`,
reached reflectively by the composite builder on the `-order META -manyOrder`
path: StateSpace died there and the decision diagram engine answered nothing
while the run looked healthy. Traced with
`Petri/test/native-trace.sh`, verified by a local rebuild, pushed as ITS-Tools
`d929bfc4` and confirmed on the cluster: the Airplane warmup on `202609080208`
is 16 examinations, 0 exceptions, `SS` answering. The config had only ever been traced over AirplaneLD PT (OneSafe,
deadlock, LTLC, UB) and COL (LTLF, CTLF, RC); it now also covers StateSpace,
Liveness, CTLCardinality, QuasiLiveness, StableMarking and
ReachabilityFireability.

## Next actions, by thread

### Native structural reductions

Inventory, pipeline and caveats: `PS_REDUCTIONS.md`; design and rules under
`Petri/src/reduction/`. Every entry point runs `reduction/Pipeline.h`'s
`prepare` (constants, initial state, drop, reduce, to a fixpoint) before any
engine; `Petri/test/reduction/check_oracle.sh MODEL EXAM [options]` is the
15 s check of one model against the deployed oracle, which is the CI archive
of 2026-09-11 with the GPPP patch, identical on the cluster.

1. Re-score `/data/ythierry/MCC26logs/local/native-reduction/*.jsonl` against
   the deployed oracle (`summarize.py`); the earlier comparisons spanned the
   broken oracle. Then rerun `mcc.py` on the pipeline binary: the answers
   before any engine and the reduced sizes changed everywhere.
2. Reduction budget: `--reductionMs` (15 s) is spent inside the run's budget
   (Erlangen bP09C09 RC: 13.7 s of reduction, nothing left to walk at 15 s).
   Share it with the engines, or bound it by the net's size.
3. STATESPACE is done for the three blocks (`reduction/Counting.h`,
   `petri64 reduce --goal STATESPACE`, `check_statespace.sh`; results in
   `PS_REDUCTIONS.md`). Open: `TRANSITIONS` after a free SCC fusion needs
   the removed internal moves per fused place and a shifted weight in
   hsc-pn (HSC_PLAN.md section 15); ITS-Tools still runs its own Java
   reductions before `hsc-pn`, the PetriSpot export is the standalone path
   until the projects integrate.
4. Dead transitions: `reduction/rules/DeadTransition.h` on
   `lp/DeadTransitions.h` runs at the head of every outer round and in the
   STATESPACE loop, `--deadMs` 3 s. Slower than libHSC's linear test on
   BugTracking by an order of magnitude (HSC_PLAN.md section 17): warm-start
   the basis between solves, add the per-place bound pre-filter, then the
   comparison and the mcc.py sweep before trusting the default budget.
5. CTL goals: a stuttering formula (no EX/AX) takes SI_CTL in Java and LTL
   here; scalar-transition removal under the LTL goal is not next-step
   preserving (Java does the same) — a theory question before enabling.
6. Review the 137 reduced timeouts and 46 limit reports of the inventory
   campaign (BlocksWorld first), 15 s per model, no cluster.
7. The blank root evidence when an until fails at its initial state (the CTL
   checker discards the child's trace) is still to repair.

### First: the CTL examinations (the two campaigns above)

Read `hsc/20260910` when it lands, then `itstools/202609080313`. The
ITS-Tools product still bundles the pre-fix `hsc-pn`: the ITS-Tools CI was
pushed (`1305abfb`, the `NATIVE_MARCH` knob) after libHSC's `HSC-Linux`
carried the fix, so the next product has it; check the pair as `docs/CI.md`
says before an ITS-Tools CTL campaign with `-hsc`.

### The LTSmin partial order soundness bug (closed on our side)

`--no-V` is back in both LTSmin runners (ITS-Tools `7f4113e0`). Without it the
reduction reports an empty product where an accepting cycle exists, and we
publish TRUE for a FALSE property: it did so once in 30 373 on LTLCardinality in
the campaign of 2026-09-07. This is LTSmin issue 169, worked around in 2019 and
commented out three days later; correct NES and NDS matrices do not make the
default visibility proviso safe. The investigation and a two command
reproduction are in `Petri/test/ltsmin-por-bug/`. What is left is upstream, and
only if someone wants it: walk the 70 step witness against the stubborn set
chosen at each of its states to find the first transition the reduction drops.

### CTL checker (`CTL_PLAN.md` section 9 has the list, 11 the design talk)

1. MCC-drivers `petrispot` tool: add CTLCardinality and CTLFireability
   (`BenchKit_head.sh` case, `SupportedExamination.txt`, README) for a
   stand-alone measurement against the oracles; a short confinement (300 s)
   for the first sweep, the checker runs to the clock on what it cannot close.
2. Saturation as a strategy whose share follows its success, on when places
   carry many tokens, never on a one-safe net (also for the reachability walk).
3. Hunts through the portfolio (strategies, quests, threads inside one
   property); the LP from a state; the evidence verifier.
4. ITS-Tools side: the Liveness pre-step spends a 30 s deadlock walk before
   the CTL loop that answers in 0.1 s (`GlobalPropertySolver`), reorder or
   shorten; the frontier handoff to `its-ctl` (section 11.3) needs an initial
   set of states as input.

### Walk engine and portfolio

1. Erlangen with 78 000 open targets after the sweep: the driver finds no
   round worth a walk and hands back 18 of 30 s; PORTFOLIO G2 (budgets from
   the clock), the sweep should go on with the time left.
2. PORTFOLIO G1: effort accounting in `DoneProperties`; the rest of the plan
   is argued from measurements the tool does not collect yet.
3. `WALK_PLAN.md` 10.10 (subgoals under a budget) is an open design; 10.11
   and 10.12 landed (tasks, coordinator). Then the families of section 9 and
   the self-configuration of section 10.
4. Watch UpperBounds and QuasiLiveness wall times: `--escalate` makes a walk
   spend its whole `--totalTime`; the knob is the `--escalate` argument in
   `PetriSpotWalker.runReachability` / `runBounds`.

### Native image (ITS-Tools)

`-march`: `build-native.sh` takes `NATIVE_MARCH` (ITS-Tools `1305abfb`,
committed, not pushed; unset keeps the CI's default). The `x86-64-v2` image
built from the deployed product is `/data/ythierry/MCC26deploy/native-v2/its-tools-native-v2`
and `cluster.lip6.fr:MCC26/flat-test/its-tools-native-v2`, tested on `small10`
and `big12` (`docs/CLUSTER.md` section 1). Open: how the deploy picks it
(`runeclipse.sh` could exec the v2 file when `/proc/cpuinfo` lacks `avx2`),
and a `small%` campaign at 6 cores. The tall probe job `1387749` is queued
until tomorrow; its log in `MCC26/MCC-drivers/probe/` gives tall's RAM per
core for the table. Open too:
`--exact-reachability-metadata` for loud misses on a sweep, and stripping the
20 signed jars at install (a 25 % start-up gain for the flat launcher). More
closed-world misses should be expected on the corpus; the recipe is
`Petri/test/native-trace.sh` then a rebuild.

### libHSC as a competitor (a dedicated session)

Exact `count`, static binaries, the other StateSpace values, deadlock, then a
StateSpace sweep (`libHSC_in_MCC.md`).

### Known and accepted

`LastZero-COL-N20` ReachabilityFireability formula 13: FALSE in the contest
and in our runs, TRUE for smpt, TAPAAL and 2025-gold; the `CPN_APPROX`
skeleton over-approximation of the coloured net, deterministic, not the
walker.

StateSpace answers three values, not four: no `TRANSITIONS`. Normal.
