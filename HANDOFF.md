# Handoff — current state, next actions

Read this first, then the design file of the thread you pick up. This file is
rewritten, never appended to: what is done leaves it (result in the README or
the design file, history in git and `docs/HISTORY.md`). State as of
2026-09-07, 21:30.

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
| the deploy tree of the harness | `/data/ythierry/MCC26deploy/MCC-drivers/` (product `202609071016` + local spotutil in `itstools/`) |
| the cluster tree | `cluster.lip6.fr:~/MCC26/MCC-drivers/`, results only there |
| result pages | `/data/ythierry/MCC26run/pages`, built by `~/git/MCC-analysis` |
| collected tables, committed | `Petri/test/mcc/csv/<campaign>/` with a README each |
| local ITS-Tools products | `/data/ythierry/itstools-local-ctl/` (this session's build, master petri64 copied in), `/data/ythierry/itstools-ci-check/` (the published product) |
| native image material | `/data/ythierry/MCC26deploy/flat-test*/`, `cluster.lip6.fr:MCC26/flat-test/` |
| development models | `bench/models/<model>/` (git-ignored; Airplane, Angiogenesis, Bridge, Erlangen, CloudDeployment-2a, CloudOpsManagement 2, 5, 10, 20, 40) |

## In flight

**Campaign 2026-09-07** (RD, QLA, LTLC, LTLF; 1800 s, 4 cores, `tall%`,
the classic product, submitted 14:37 CEST). At 18:35: RD complete and
collected (1953 logs, 0 wrong, variance against 06c), QLA 1681 of 1953, LTLC
1222, LTLF not started, 111 jobs running, 2686 waiting. Collect QLA, LTLC and
LTLF when they end: `rsync -rz --exclude='*.stderr' cluster.lip6.fr:MCC26/MCC-drivers/{QLA,LTLC,LTLF} /data/ythierry/MCC26run/2026-09-07/`
(never `--delete`), then `mcclogs2csv.py LTLC LTLF -o csv`,
`totallogs2csv.py QLA -o csv`, `report.py csv --baseline
Petri/test/mcc/csv/2026-09-06-baseline` for LTL (the question: did the Effort
glean bring LTL back to the contest-like ITS-Tools? Airplane says yes) and
`csv/2026-09-06c` for QLA; point `ITS-Tools latest` in
`~/git/MCC-analysis/campaign/example.json` at the folder, rebuild the pages,
archive, remove from the cluster. Expect about 5 Eclipse fatals per
examination; check the `failure` column before believing a one second miss.

**CI.** PetriSpot `Inv-Linux` carries the CTL checker (deployed from 16cb816);
the published ITS-Tools product `202609071619` bundles that binary (checksum
verified) and the CTL plug. ITS-Tools 715a2036 (a message fix) is pushed and
its product is due; nothing on the cluster uses either yet (the campaign
holds `itstools/`).

**Cluster bench jobs** 1365883 (tall) and 1365884 (small) time the three
launchers (Eclipse, flat, native image), queued behind the campaign.

## Next actions, by thread

### First: the CTL examinations on the cluster, once campaign 2026-09-07 has drained

The explicit CTL checker (`CTL_PLAN.md`) is built, plugged into ITS-Tools
beside `its-ctl` and published by both CIs; it has never run at scale. The
campaign holds `itstools/` on the cluster, so wait for `cluster_status.sh` to
show LTLF complete and the queue empty, collect and archive (`docs/CLUSTER.md`
sections 4 and 5), then:

1. Deploy the published product (`docs/CLUSTER.md` section 1): fresh
   `install_itstools.sh` in the deploy tree, check that the product's bundled
   `petri64 -h` lists `--ctlSteps` and that its stamp is `202609071619` or
   later (`docs/CI.md`), rsync the `itstools/` subtree.
2. Warm up on `oracle/AirplaneLD-PT-0010-*.out`; in the CTLC and CTLF logs
   look for `CTL check beside the decision diagrams` and verdicts tagged
   `CTL_WALK` (locally: CTLF 11 of 16, CTLC 5 of 16 by the checker, all 32
   right, `Petri/test/logs/its-airplane-ctl{c,f}.log`).
3. Submit CTLC and CTLF, 1800 s, 4 cores, `tall%`, one after the other
   (`submit-<date>.sh` pattern, `Petri/test/mcc/submit-2026-09-07.sh`), and
   Liveness after them if time allows: it takes the same path
   (`GlobalPropertySolver` states it as `AG EF fireable` per transition) and
   the campaign of 09-05 had it at 48 % of its wall time in the diagram tail.
4. Collect with `collect.sh <date> CTLC CTLF`; there is no complete CTL
   baseline of ours (the 09-06 CTL runs were deleted), so read the report
   against the field (`report.py --raw` for the tool board, the pages against
   `ITS-Tools 2026`, `Tapaal 2026`): the questions are how many formulas the
   checker answers (`grep -c CTL_WALK` over the logs, and how many of them
   only the checker had before the diagrams), any wrong verdict (a soundness
   bug; the oracle script found two locally, none remain), and whether the
   companion costs anything on the instances the diagrams solved alone.
5. Optionally, the stand-alone measurement: the `petrispot` tool of
   MCC-drivers gets CTLC and CTLF (below), 300 s confinement for a first
   sweep.

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

After the campaign: `install_itstools.sh` in a fresh `itstools/` on the
cluster (or the image dropped beside `its-tools`), the Airplane warmup
through `run_oar.sh`, read the logs for closed-world misses (named
exceptions), then one examination. If `small` refuses the image, `-march` in
`build-native.sh`. Open: `--exact-reachability-metadata` for loud misses on
a sweep; stripping the 20 signed jars at install (a 25 % start-up gain for
the flat launcher).

### libHSC as a competitor (a dedicated session)

Exact `count`, static binaries, the other StateSpace values, deadlock, then a
StateSpace sweep (`libHSC_in_MCC.md`).

### Known and accepted

`LastZero-COL-N20` ReachabilityFireability formula 13: FALSE in the contest
and in our runs, TRUE for smpt, TAPAAL and 2025-gold; the `CPN_APPROX`
skeleton over-approximation of the coloured net, deterministic, not the
walker.
