# Handoff — current state, next actions

Read this first, then the design file of the thread you pick up. This file is
rewritten, never appended to: what is done leaves it (result in the README or
the design file, history in git and `docs/HISTORY.md`). State as of
2026-09-08, 05:35.

## Orientation

| thread | read | code |
| --- | --- | --- |
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
| campaign logs, collected | `/data/ythierry/MCC26logs/<tool>/<build>/<EXAM>/`, `csv/` beside them; `_local/<name>/` a local experiment's files |
| the deploy tree of the harness | `/data/ythierry/MCC26deploy/MCC-drivers/` (product `202609080208`) |
| the cluster tree | `cluster.lip6.fr:~/MCC26/MCC-drivers/`, results only there |
| result pages | `/data/ythierry/MCC26logs/_shared/pages`, built by `~/git/MCC-analysis` |
| the published product, unzipped for a check | `/data/ythierry/MCC26deploy/products/itstools-ci-check/` |
| native image material | `~/git/ITStools/ITS-commandline/native/`, `cluster.lip6.fr:MCC26/flat-test/`, local test image `/data/ythierry/MCC26deploy/its-tools-native-test` |
| GraalVM for a local image build | `/data/ythierry/graal/graalvm-jdk-25.0.4+7.1` |
| development models | `bench/models/<model>/` (git-ignored) |

## In flight

**Campaign 202609080313, the first CTL examinations at scale.** Running:
launched 2026-09-08 05:33 by `Petri/test/mcc/submit-2026-09-08.sh`, CTLC then
CTLF then Liveness, 1800 s, 4 cores, `tall%`, on the native image of the product
`202609080313`, 5859 jobs, about 13 hours. That product carries both fixes of
the night, `--no-V` and the coloured `Sort[]`. Collect into a folder named for
it:

```
bash Petri/test/mcc/collect.sh 202609080313 CTLC CTLF L
```

then a new set in `~/git/MCC-analysis/campaign/example.json`, rebuild the pages,
archive, free the cluster (`docs/CLUSTER.md` sections 4, 5). There is no CTL
baseline of ours, so read it against the field (`report.py --raw`, the pages
against `ITS-Tools 2026` and `Tapaal 2026`): how many formulas the checker
answers, how many only it had before the diagrams, any wrong verdict, and
whether the companion costs anything where the diagrams won alone. Liveness
takes the CTL path too.

A first attempt went out at 04:47 on `202609080208` and was cancelled an hour
later, all 3067 jobs deleted and the folders cleared: every coloured model
declaring a product sort died on `symmetricnet.terms.Sort[]`, answering nothing.
Should a handful fail this time, nothing needs redoing wholesale --
`Petri/test/mcc/resubmit.sh <EXAM>` submits only what has no log, did not reach
the teamcity suite close, or carries a closed world miss.

The warmup before this launch covered AirplaneLD-PT-0010 **and BART-COL-002**,
all 16 examinations, everything answering, no exception; CTLC 5 and CTLF 24
verdicts tagged `CTL_WALK`. Keep the coloured instance in the warmup: Airplane
alone passed the aborted campaign's warmup. Note also that the published image
was byte for byte the same *size* as the broken one, so a native image is
checked by running it, never by its size.

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
