# Handoff — current state, next actions

Read this first, then the design file of the thread you pick up. This file is
rewritten, never appended to: what is done leaves it (result in the README or
the design file, history in git and `docs/HISTORY.md`). State as of
2026-09-10, 14:30.

## Orientation

| thread | read | code |
| --- | --- | --- |
| the walk engine (reachability) | `WALK_PLAN.md` sections 9, 10 | `Petri/src/walk/`, `Petri/src/sched/` |
| the solving loop, walker budgets | `PORTFOLIO.md`, `INTEROP.md` | `cli/WalkDriver.h`, ITS-Tools `PetriSpotWalker` |
| the CTL checker | `CTL_PLAN.md` sections 9, 11 | `Petri/src/ctl/`, `cli/CtlDriver.h` |
| the state equation, LP hints | `Petri/src/lp/algorithm.md` | `Petri/src/lp/` |
| the campaigns, the harness | `BENCH.md`, `Petri/test/mcc/README.md`, `TOTAL_QUERIES.md` | `Petri/test/mcc/` |
| libHSC as a competitor and companion | `HSC_PLAN.md`, `HSC_EXPERIMENTS.md`, `libHSC_in_MCC.md`; state in libHSC `handoff_mcc.md` | `~/git/MCC-drivers/hsc/`, ITS-Tools `hsc/`, `interop/` |
| the libHSC order sweep | libHSC `experiments/order/SWEEP.md` §9, `handoff_order.md` | libHSC `include/hsc/order/` |
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
| collected campaign logs | `/data/ythierry/MCC26logs/<tool>/<build>/`, `csv.<build>/` beside them |
| the deploy tree of the harness | `/data/ythierry/MCC26deploy/MCC-drivers/` |
| the cluster tree | `cluster.lip6.fr:~/MCC26/MCC-drivers/`, results only there |
| the order sweep on the cluster | `cluster.lip6.fr:~/MCC26/hsc-sweep/`, reaped to `MCC26logs/hsc/sweep1/` |
| result pages | `/data/ythierry/MCC26logs/_shared/pages`, built by `~/git/MCC-analysis` |
| collected tables, committed | `Petri/test/mcc/csv/<campaign>/` with a README each |
| local ITS-Tools products | `/data/ythierry/itstools-ci-check/` (the published product) |
| native image material | `~/git/ITStools/ITS-commandline/native/`, `cluster.lip6.fr:MCC26/flat-test/`, local test image `/data/ythierry/MCC26deploy/its-tools-native-test` |
| GraalVM for a local image build | `/data/ythierry/graal/graalvm-jdk-25.0.4+7.1` |
| development models | `bench/models/<model>/` (git-ignored) |

## In flight

**Nothing runs on the cluster and our queue is empty.** The 2651 waiting
jobs of the libHSC order sweep were cancelled on 2026-09-10: every `tall`
node is held by other users' whole-node jobs and the scheduler's estimate
for our next job was four days out. `tall%` is still full; `small%` (24
nodes, 576 cores) and `big%` are idle, and `docs/CLUSTER.md` §1a says what
`small%` is and what runs there.

**The campaign of 2026-09-10 is ready and unsubmitted.** Everything is
deployed, pushed and warmup-tested locally; the plan, the legs and what was
tested are in `Petri/test/mcc/campaign-2026-09-10.md`, the script is
`Petri/test/mcc/submit-2026-09-10.sh`. Legs 2 and 3 (`hsc` CTLC/CTLF at
1800 s, `hscapprox` RC/RF at 300 s) can go to `small%` today: `hsc-pn` runs
there unmodified. Leg 1 (`itstools` RC/RF) needs `tall%` until the native
image is rebuilt at `-march=x86-64-v2`.

## Next actions, by thread

### First: a verdict off a partial reachable set is unsound (libHSC)

The order sweep found it before it found anything about shapes
(libHSC `experiments/order/SWEEP.md` §9.1, `handoff_ctl.md` item 0). The
CTL checker answers from a truncated `R`: on `FMS-PT-10000` under `louvain`
the run prints `R partial`, `partial=1`, `reach_states=8764`, then answers
`FALSE` on a `TRUE` formula. Across the sweep, 545 of the 559 wrong CTLC
verdicts come from a run whose `R` was partial, and **47 % of every answer
given on a partial `R` was wrong**. It was invisible until now because
every earlier CTL measurement closed `R`. This invalidates the "0 wrong" of
every budgeted libHSC CTL run, and it must be fixed before leg 2 of the
campaign means anything.

Separately, 14 wrong verdicts sit on a *complete* `R` (§9.2) — a shape
changed an answer on a full state space, which nothing excuses. Six
instances, listed there.

### The CTL examinations at scale (ITS-Tools) — collected and read

Campaign `202609080313`, 1954 instances each at 1800 s, in
`/data/ythierry/MCC26logs/itstools/202609080313/`:

| | answered | ok | wrong | bonus |
| --- | ---: | ---: | ---: | ---: |
| CTLCardinality | 26069 | 25051 | **1** | 1017 |
| CTLFireability | 23682 | 22579 | **2** | 1101 |
| Liveness | 1814 | 1811 | 0 | 3 |

The three wrong verdicts are the open item, none backed by another tool:

| formula | oracle | ours | run s |
| --- | --- | --- | ---: |
| `ShieldIIPs-PT-002A-CTLCardinality-2024-07` | FALSE | TRUE | 307 |
| `ClientsAndServers-PT-N0020P1-CTLFireability-2024-09` | FALSE | TRUE | 1801 |
| `FileSystem-COL-N02I10B10-CTLFireability-2024-05` | TRUE | FALSE | 767 |

Two are TRUE for a FALSE property and one the converse, so no single
direction of unsoundness explains them; `FileSystem-COL` is coloured, which
puts the skeleton approximation on the list of suspects for that one alone.
Three wrong in 49 565 answers. The tables are committed at
`Petri/test/mcc/csv/202609080313/`, its README the read.

### The native image (ITS-Tools)

`build-native.sh` needs `-march=x86-64-v2`, which covers `tall%` and
`small%` in one image; `-march=compatibility` is the fallback that cannot
fail. The failure is recorded at `cluster.lip6.fr:~/MCC26/flat-test/run-native.log`
and quoted in `docs/CLUSTER.md` §1. This is the one flag between us and
half the cluster. Open too: `--exact-reachability-metadata` for loud misses
on a sweep, and stripping the 20 signed jars at install (a 25 % start-up
gain for the flat launcher). More closed-world misses should be expected on
the corpus; the recipe is `Petri/test/native-trace.sh` then a rebuild.

### The LTSmin partial order soundness bug (closed on our side)

`--no-V` is back in both LTSmin runners (ITS-Tools `7f4113e0`). Without it
the reduction reports an empty product where an accepting cycle exists, and
we publish TRUE for a FALSE property: it did so once in 30 373 on
LTLCardinality in the campaign of 2026-09-07. This is LTSmin issue 169,
worked around in 2019 and commented out three days later; correct NES and
NDS matrices do not make the default visibility proviso safe. The
investigation and a two command reproduction are in
`Petri/test/ltsmin-por-bug/`. What is left is upstream, and only if someone
wants it: walk the 70 step witness against the stubborn set chosen at each
of its states to find the first transition the reduction drops.

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

### libHSC as a competitor (a dedicated session)

Exact `count`, static binaries, the other StateSpace values, deadlock, then a
StateSpace sweep (`libHSC_in_MCC.md`).

### Known and accepted

`LastZero-COL-N20` ReachabilityFireability formula 13: FALSE in the contest
and in our runs, TRUE for smpt, TAPAAL and 2025-gold; the `CPN_APPROX`
skeleton over-approximation of the coloured net, deterministic, not the
walker.

StateSpace answers three values, not four: no `TRANSITIONS`. Normal.
