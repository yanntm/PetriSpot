# Handoff — current state, next actions

Read this first, then the design file of the thread you pick up. This file is
rewritten, never appended to: what is done leaves it (result in the README or
the design file, history in git and `docs/HISTORY.md`). State as of
2026-09-08, 04:10.

## Orientation

| thread | read | code |
| --- | --- | --- |
| the walk engine (reachability) | `WALK_PLAN.md` sections 9, 10 | `Petri/src/walk/` |
| the solving loop, walker budgets | `PORTFOLIO.md`, `INTEROP.md` | `cli/WalkDriver.h`, ITS-Tools `PetriSpotWalker` |
| the CTL checker | `CTL_PLAN.md` sections 9, 11 | `Petri/src/ctl/`, `cli/CtlDriver.h` |
| the state equation, LP hints | `Petri/src/lp/algorithm.md` | `Petri/src/lp/` |
| the campaigns, the harness | `BENCH.md`, `Petri/test/mcc/README.md`, `TOTAL_QUERIES.md` | `Petri/test/mcc/` |
| libHSC as a competitor | `libHSC_in_MCC.md` | `~/git/MCC-drivers/hsc/` |
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
| campaign logs, collected | `/data/ythierry/MCC26run/<date>/<EXAM>/`, `csv/` beside them |
| archived campaigns | `/data/ythierry/MCC26archive/<date>/` |
| the deploy tree of the harness | `/data/ythierry/MCC26deploy/MCC-drivers/` (product `202609071637`) |
| the cluster tree | `cluster.lip6.fr:~/MCC26/MCC-drivers/`, results only there |
| result pages | `/data/ythierry/MCC26run/pages`, built by `~/git/MCC-analysis` |
| collected tables, committed | `Petri/test/mcc/csv/<campaign>/` with a README each |
| local ITS-Tools products | `/data/ythierry/itstools-ci-check/` (the published product) |
| native image material | `~/git/ITStools/ITS-commandline/native/`, `cluster.lip6.fr:MCC26/flat-test/`, local test image `/data/ythierry/MCC26deploy/its-tools-native-test` |
| GraalVM for a local image build | `/data/ythierry/graal/graalvm-jdk-25.0.4+7.1` |
| development models | `bench/models/<model>/` (git-ignored) |

## In flight

**The CTL examinations, blocked on one CI run.** The cluster is free (25 GB,
inputs only; campaign 2026-09-07 archived and its result folders removed). The
plan is CTLC, CTLF then Liveness at 1800 s, 4 cores, `tall%`, on the **native
image**, by `Petri/test/mcc/submit-2026-09-08.sh` (5859 jobs, about 13 hours).

Waiting on the ITS-Tools CI for `d929bfc4` (the reachability metadata fix
below). When it is green: reinstall the product in the deploy tree, check the
stamp is later than `202609071637`, rsync `itstools/`, re-run the Airplane
warmup and check that `SS` answers its three values, then submit.

**The native image is now the launcher.** `runeclipse.sh` execs
`its-tools-native` whenever the file is present, so the deploy decides it;
`docs/CLUSTER.md` section 1 says how to tell and how to choose. The image is
AVX2, hence `tall%` only. Its closed world was missing
`fr.lip6.move.gal.InstanceDecl[]` and `fr.lip6.move.gal.Synchronization[]`,
reached reflectively by the composite builder on the `-order META -manyOrder`
path: StateSpace died there and the decision diagram engine answered nothing
while the run looked healthy. Traced with
`Petri/test/native-trace.sh`, verified by a local rebuild, pushed as ITS-Tools
`d929bfc4`. The config had only ever been traced over AirplaneLD PT (OneSafe,
deadlock, LTLC, UB) and COL (LTLF, CTLF, RC); it now also covers StateSpace,
Liveness, CTLCardinality, QuasiLiveness, StableMarking and
ReachabilityFireability.

## Next actions, by thread

### First: the CTL examinations (the campaign above)

The explicit CTL checker (`CTL_PLAN.md`) has never run at scale. Through the
native image the Airplane warmup gives CTLC 5 of 16 and CTLF 11 of 16 verdicts
tagged `CTL_WALK`, all right, 0 exceptions -- the local numbers exactly.

Collect with `collect.sh 202609071637 CTLC CTLF L` (the run folder is named
after the product stamp). There is no complete CTL baseline of ours, so read
the report against the field (`report.py --raw`, the pages against
`ITS-Tools 2026`, `Tapaal 2026`): how many formulas the checker answers, how
many of those only it had before the diagrams, any wrong verdict (a soundness
bug), and whether the companion costs anything where the diagrams won alone.

### The LTSmin partial order soundness bug

`StigmergyCommit-PT-02b-LTLCardinality-03` is answered TRUE where the field has
FALSE, by `PARTIAL_ORDER EXPLICIT LTSMIN SAT_SMT`, and it was answered FALSE on
2026-09-06 by the same technique on the same net. The knowledge step is what
changed (four factoids that reduced the automaton 3 states/4 edges to 2/3 on
09-06 reduce nothing on 09-07). Weaker knowledge must cost an answer, not
produce a wrong one. Details in `Petri/test/mcc/csv/2026-09-07/README.md`.

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

`-march`: the image refuses to start on `small%` (no AVX2) and answers nothing
rather than failing loudly, so half the cluster is out of reach. Open too:
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
