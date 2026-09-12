# Handoff — current state, next actions

Read this first, then the design file of the thread you pick up. This file is
rewritten, never appended to: what is done leaves it (result in the README or
the design file, history in git and `docs/HISTORY.md`). State as of
2026-09-12, 01:00 CEST.

## Orientation

| thread | read | code |
| --- | --- | --- |
| standing alone: examinations, procedures, the loop | `AUTONOMY.md` (design, to comment) | none yet |
| native structural reductions | `PS_REDUCTIONS.md`, `Petri/src/reduction/README.md`, `Petri/src/reduction/algorithm.md` | `Petri/src/reduction/`, `Petri/test/reduction/` |
| the walk engine (reachability) | `WALK_PLAN.md` sections 9, 10 | `Petri/src/walk/` |
| the solving loop, walker budgets | `PORTFOLIO.md`, `INTEROP.md` | `cli/WalkDriver.h`, ITS-Tools `PetriSpotWalker` |
| the CTL checker | `CTL_PLAN.md` sections 9, 11 | `Petri/src/ctl/`, `cli/CtlDriver.h` |
| the state equation, LP hints | `Petri/src/lp/algorithm.md` | `Petri/src/lp/` |
| structural inequalities and projection bounds | `handoff_inequalities.md`, `INEQUALITIES.md` | `Petri/src/invariants/`, libHSC `tools/pn_approx*.hh` |
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
| the deploy tree of the harness | `/data/ythierry/MCC26deploy/MCC-drivers/` (HSC/PetriSpot updated; ITS product separately staged) |
| the cluster tree | `cluster.lip6.fr:~/MCC26/MCC-drivers/`, results only there |
| result pages | `/data/ythierry/MCC26logs/web/campaign`, built by `~/git/MCC-analysis`; `web/order-sweep` the libHSC sweep's |
| native image material | `~/git/ITStools/ITS-commandline/native/`, `cluster.lip6.fr:MCC26/flat-test/`, the x86-64-v2 image `/data/ythierry/MCC26deploy/native-v2/` |
| GraalVM for a local image build | `/data/ythierry/graal/graalvm-jdk-25.0.4+7.1` |
| development models | `bench/models/<model>/` (git-ignored) |

## In flight

**Reduced PT campaigns: detached submission running.** Controller:
`cluster.lip6.fr:~/MCC26/MCC-drivers/submit-reduce-20260912.sh full`,
log `submit-reduce20260912.log` beside it. Do not rewrite running scripts
or deployed tools. One serial controller; no parallel oarsub. Queue checks
only between full sweep batches, every five minutes when capacity is needed;
no queue drain between batches. Up to 5000 queued jobs is acceptable.

All 1681 PT models, no COL oracles. StateSpace: 900 s, small%, six cores,
17-minute OAR walltime. Then the original shape sweep, SS/CTLC/CTLF,
17 heuristics per bundled job, proven duplicate net+shape pairs reused;
300 s external / 270 s internal, original 94-minute bundle walltime.
Sweep deployment: `~/MCC26/hsc-sweep-reduce-20260912/`, results under
`results/reduce20260912/`. StateSpace logs: `SS.reduce20260912/` in the harness.
Total planned jobs: 1681 StateSpace + 5043 sweep bundles.

Warmups reviewed: three StateSpace jobs, six sweep jobs. All 102 sweep entries
accounted for, zero wrong answers or overruns. StateSpace returned 11 correct
values; DoubleExponent omitted TRANSITIONS, the accepted reduction limitation.
Local evidence: `/data/ythierry/MCC26logs/hsc/reduce20260912/`.
Use the sweep audit with the explicit model and heuristic rosters when collecting.

**Follow-up after submission:** confirm the controller eventually reaches
SUBMISSION DONE; collect campaign results locally, merge sweep rows, rebuild
pages. No broad cluster-side analysis. RC/RF build 202609101624 is already
collected (3906 logs) and its CSV/report rebuilt under
`/data/ythierry/MCC26logs/itstools/202609101624/`; its page configuration,
page rebuild and log-index updates remain. Do not collect it again unnecessarily.
Older HSC campaign collection completeness was not rechecked this session.

Counting contract: `Petri/src/io/PNET.md` and libHSC's tool documentation.
PCONST removes constant free components from the DD and preserves exact GMP
binomial factors; live PCOEF stays in DD counting. Transition reconstruction
for these reductions remains open. Current paired CI binaries are deployed;
source changes pushed in PetriSpot 6749234, libHSC bbc823b, harness 976ff018,
sweep tooling 8532db4. The sweep records raw paths, DD nodes, exact weighted
markings, constant factor, epochs and partial completion.

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
3. `hsc-pn --reduce` (libHSC, HSC_PLAN.md section 17) runs these
   reductions in memory before the portfolio; the reduced StateSpace
   campaign is in flight above. Counting metadata is documented in PNET.md (`reduction/Counting.h`,
   `petri64 reduce --goal STATESPACE`, `check_statespace.sh`; results in
   `PS_REDUCTIONS.md`). Open: `TRANSITIONS` after a free SCC fusion needs
   the removed internal moves per fused place and a shifted weight in
   hsc-pn (HSC_PLAN.md section 15); ITS-Tools still runs its own Java
   reductions before `hsc-pn`, the PetriSpot export is the standalone path
   until the projects integrate.
4. Dead transitions: `reduction/rules/DeadTransition.h` on
   `lp/DeadTransitions.h` runs at the head of every outer round and in the
   STATESPACE loop, `--deadMs` 3 s, places first then transitions. On
   BugTracking libHSC's invariant-set test finds more in 0.2 s than the LP
   in 30 s, an integrality gap not a budget one (HSC_PLAN.md section 17):
   hsc-pn computes the invariants and tests on the set; PetriSpot keeps
   the LP under its budget, to be swept with mcc.py before trusting the
   default on every examination.
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

Local Maven and x86-64-v2 native build completed, staged at
`/data/ythierry/MCC26deploy/products/itstools-202609112344/`.
ITS deployment was deprioritized for tonight's HSC-only campaigns. That product
predates the newest PCONST binaries: refresh its bundled HSC/PetriSpot before
using it for this counting contract. Do not assume it replaced the cluster image.

Side note: signatures were stripped from 20 product JARs to get past a native
build certificate-metadata failure. Harmless workaround accepted for this run;
the signed-JAR issue and whether signatures were verified at startup remain
**uninvestigated**. Do not present this as a diagnosed startup issue.

### libHSC as a competitor (a dedicated session)

Exact `count`, static binaries, the other StateSpace values, deadlock, then a
StateSpace sweep (`libHSC_in_MCC.md`).

### Known and accepted

`LastZero-COL-N20` ReachabilityFireability formula 13: FALSE in the contest
and in our runs, TRUE for smpt, TAPAAL and 2025-gold; the `CPN_APPROX`
skeleton over-approximation of the coloured net, deterministic, not the
walker.

Reduced HSC StateSpace may omit `TRANSITIONS` when reconstruction is unsupported;
it reports that value when preserved. Missing it alone is not a wrong answer.
