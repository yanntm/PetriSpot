# Handoff — MCC 2026 campaign, the portfolio design, the first fixes, the total examinations

State of the work as of 2026-09-06, 11:30. Read `PORTFOLIO.md` first if you are
picking up the design, `TOTAL_QUERIES.md` for the total examinations; read
this for where everything lives and what is in flight.

## 2026-09-07, 17:30: closing the session's threads

**Campaign 2026-09-07, 17:45: RD collected.** `/data/ythierry/MCC26run/2026-09-07/RD`
(1953 logs, `csv/` beside it): 1897 answered as in 06c, 0 wrong, 3 gained
(DLCflexbar-PT-8b, PGCD-COL-D02N100, RERS2020-PT-pb104), 3 lost
(DatabaseWithMutex-PT-20, HypercubeGrid-PT-C5K3P3B15,
HypertorusGrid-PT-d5k3p2b10), the same 5 Eclipse fatals and 2 overlarge
markings: variance. The pages (`/data/ythierry/MCC26run/pages`, built 16:37)
carry it as the set `ITS-Tools 2026-09-07`; the previous `latest` is now named
`ITS-Tools 2026-09-06c` (MCC-analysis c3aae82). QLA, LTLC, LTLF follow the
recipe below when they end. The deploy clone of MCC-drivers is repaired and at
origin (b6440228): the 5 050 untracked collisions were byte-identical copies.

**Campaign 2026-09-07 (RD, QLA, LTLC, LTLF), running.** Submission complete
at 14:37 CEST (`submit-2026-09-07.log`); at 17:00 RD had all 1953 logs, QLA
1195, LTLC and LTLF none yet, 4505 jobs queued, about 110 running. To
analyse it: `rsync -rz --exclude='*.stderr' cluster.lip6.fr:MCC26/MCC-drivers/{RD,QLA,LTLC,LTLF} /data/ythierry/MCC26run/2026-09-07/`
(without `--delete`; the stderr carry the harness `time -p` lines if wanted),
then `mcclogs2csv.py RD LTLC LTLF -o csv` and `totallogs2csv.py QLA -o csv`,
`report.py csv --baseline ~/git/PetriSpot/Petri/test/mcc/csv/2026-09-06-baseline`
for LTLC and LTLF (the contest-like ITS-Tools; the question is whether the
Effort glean brought LTL back to it: Airplane says yes), `csv/2026-09-06c`
for RD and QLA (the walker of 6a5d288 against 3f29014 with the step cap);
point `ITS-Tools latest` in `~/git/MCC-analysis/campaign/example.json` at the
new folder and rebuild the pages; archive whole to
`/data/ythierry/MCC26archive/2026-09-07/` and remove from the cluster as done
for 06c. Watch the first-wave Eclipse fatals (5 in RD, the usual) and any run
whose log names a class: none expected, the product is the classic one.

**Native image, state.** Published by the CI (`its-tools-native`, 80 MB,
14:20 UTC, the beside-the-executable lookup), driver support pushed
(ITS-Tools-MCC 3839363), verified locally through `run_test.pl` (Airplane
LTLC 16/16 in 10.7 s, RD 0.05 s). Not yet run on the cluster: the campaign
holds `itstools/`; after it, `install_itstools.sh` in a fresh clone there
(or the image dropped into `itstools/itstools/` beside `its-tools`), the
Airplane warmup, then one examination. Bench jobs 1365883 (tall) and 1365884
(small) still queued; if `small` refuses the image, `-march` in
`build-native.sh`. Open items on that path: `--exact-reachability-metadata`
for loud closed-world misses on a sweep, stripping the 20 signed jars at
install (a 25 % start-up gain for the flat launcher, irrelevant to the image).

**libHSC as a competitor.** Explored in `libHSC_in_MCC.md` (this repository)
with the measurements; a portfolio driver prototype `MCC-drivers/hsc/` is
committed and pushed (b6440228) and works through `run_test.pl`. Work for a
dedicated session: exact `count`, static binaries, the other StateSpace
values, deadlock, then a StateSpace sweep. The deploy clone
`/data/ythierry/MCC26deploy/MCC-drivers` is behind origin: untracked local
oracle files collide with ones upstream now tracks; move them aside and pull.

## 2026-09-07, 16:50: the image published, the site rewritten

ITS-Tools CI: the native step and the Node 24 update (`checkout@v7`,
`setup-java@v6`, the runner's own Maven; c6f632b0, 043bd1ca) both green,
`https://lip6.github.io/ITSTools/its-tools-native` published (80 MB).
54c9dbac makes the image find `plugins/` beside its own executable (the code
source of a native image is the executable; a file that is not a jar means
that) so no property is needed, and rewrites `website/index.html` and the
`README.md`: Java 21, the `-pnfolder`/`-examination` and `-i model.gal`
command lines, the flat script, the native image as an experimental build
made for the contest driver with how to report a closed-world miss and the
fallback to `its-tools`. Verified locally before the push: image beside
`plugins/`, no property, from another directory: LTLC 16/16, COL LTLF 16/16,
deadlock, all binaries found. The cluster's `itstools/` is not touched while
the campaign runs; the local deploy tree gets the published image through the
updated driver (`git pull` in `MCC-drivers/itstools`, `wget` of the image)
for a `run_test.pl` check, then the cluster trial after the campaign.

## 2026-09-07, 16:10: the native image wired into the CI and the MCC driver

ITS-Tools (commit after 3aec3d8c): `main()` in both application classes runs
the tool on a thread with a 128 MB stack (`Application.MAIN_STACK`, the
ini's `-Xss128m`; a native executable's `main` would otherwise sit on the OS
stack, 64 MB under the driver's `ulimit -s 65536`); `build.yml` adds
`graalvm/setup-graalvm` (Java 25, Oracle GraalVM) and runs
`ITS-commandline/native/build-native.sh` on the Linux product (`NATIVE_XMX`
10g for the 16 GB runner), publishing `its-tools-native` on gh-pages beside
the product zip; the travis-era `ITS-commandline/{runeclipse.sh,
install_eclipse.sh,.travis.yml}` are removed (the driver has its own
`runeclipse.sh`). Verified locally before the push: product rebuilt, image
rebuilt from it in 36 s, OneSafe 18 ms, LTLC and COL LTLF 16/16, no error.
ITS-Tools-MCC 3839363 (`~/git/MCC-drivers/itstools`, pushed):
`install_itstools.sh` fetches `its-tools-native` after the zip (dropped when
the URL fails), `runeclipse.sh` execs it with
`-Dfr.lip6.binaries.root=$BINDIR/itstools/plugins` when present, the
Eclipse launcher otherwise. Meant for a sweep: if the image holds, one of the
two products suffices. Bench jobs on the cluster: 1365883 (tall) and 1365884
(small, the old hardware that may lack AVX2; `-march` in build-native.sh is
the knob), both `~/MCC26/flat-test/flatbench.sh`, queued behind the campaign.
Next: when the CI has published the image, `install_itstools.sh` in a fresh
`MCC-drivers/itstools/` on hydrogen, the Airplane warmup through
`run_oar.sh`, read the logs for closed-world misses (named exceptions), trace
and rebuild if any, then a real examination.

## 2026-09-07, 15:30: ITS-Tools as a native executable (GraalVM), it works

Oracle GraalVM 25.0.4 in `/data/ythierry/graal/` (no Fedora package; the
`zlib-ng-compat-devel` already installed provides zlib). The tracing agent
over seven examinations of AirplaneLD PT-0010 and COL-0010 recorded a small
closed world (182 reflective types, 35 ours, 30 resources); `native-image`
over the flat class path built a 60 MB executable in 36 s (77 MB with G1).
Ten examinations against the flat JVM launcher, **identical verdicts every
time, the two untraced examinations (CTLC, RF) included**: OneSafe 10 ms
(flat 240 ms, Eclipse 340 ms warm, 1.2 s cold); LTLC 8.4 s vs 9.3 s; CTLC
4.5 vs 5.1; COL LTLF 1.1 vs 2.2; RF and UB equal (75 and 105 s in the
external engines). No JIT-less slowdown seen; the missing warm-up shows as a
gain. The G1 build (`--gc=G1 -R:MaxHeapSize=16g -R:MinHeapSize=40m
-R:StackSize=128m`, the ini's values as run-time defaults) runs the same.
ITS-Tools `08a3de1b` gives the binaries plugins
`-Dfr.lip6.binaries.root=<product>/plugins` (an image has no class folder);
`ITS-commandline/native/` (commit after it) holds `build-native.sh`,
`trace-agent.sh`, the recorded `config/reachability-metadata.json` and a
README. A closed-world miss at run time is a `ClassNotFoundException`,
`NoSuchMethodException` or `MissingReflectionRegistrationError` naming the
element (a resource: null, a few frames later); trace that run, rebuild.
The image `its-tools-native` (G1) sits beside the flat product in
`cluster.lip6.fr:MCC26/flat-test/`, and job 1365883 (`flatbench.sh`, still
queued behind the campaign) now times the three launchers on a tall node.
A cluster trial of the image on real examinations is a variant of
`runeclipse.sh` (ITS-Tools-MCC) calling the executable with the binaries
root; not done. Local artefacts: `/data/ythierry/graal/{build.sh,build.log,
build-g1.log,runs/}`, products `/data/ythierry/MCC26deploy/flat-test{2,3}/`.

## 2026-09-07, 14:50: ITS-Tools on a flat classpath, a prototype

**ITS-Tools f9e46d4e.** `its-tools-flat.sh`, shipped beside `its-tools` in the
Linux product (a root file of `fr.lip6.move.gal.itscl.feature`): `java -cp`
over the product's `plugins/` (every jar, every unpacked bundle folder, the
jars a bundle nests under `lib/` unpacked into a private temp folder), the
JVM flags of `its-tools.ini`, main class
`fr.lip6.move.gal.itscl.application.Application`. No Equinox, no launcher, no
configuration area, no extension registry: the "Cannot open display" fatals
of the first wave cannot happen. Code: `main(String[])` in both application
classes, the `IApplication.start` delegating to a `run(args)`; the six
binaries plugins locate `bin/` from the folder their class was loaded from
when `getDefault()` is null (no framework), the Eclipse path unchanged.
EMF and Xtext need nothing: `setStandalone(true)` was already on for the CLI
and the generated packages self-register. Verified against the Eclipse
launcher on AirplaneLD-PT-0010 LTLC, deadlock and OneSafe, AirplaneLD-COL-0010
RC (the PNML framework, unfolding, skeleton) and LTLF (spotutil): same
verdicts, same times, no error. The stashed BND experiment
(`stash@{0}`, a bnd-export of a bndrun) was the other road: one executable
jar still launching Equinox and `eclipse.application`; it stalled on bnd's
FileRepo not reading Tycho's flat `plugins/` layout, and is not needed here.

**Startup, measured.** Locally, OneSafe on AirplaneLD-PT-0010: Eclipse
launcher 340-400 ms wall for 41-44 ms of `Total runtime`; flat 310-600 ms
wall for 97-203 ms internal (a 75-entry flat classpath makes every first
class load a linear jar search, where Equinox indexes packages), a bare JVM
10 ms. Through `run_test.pl` locally the gap wall minus internal is 0.77 s.
**On the cluster** (the 260 Airplane warmup runs, `time -p` in the stderr
against our `Total runtime`): gap median 1.30 s, p10 1.02, p90 1.74, never
under 0.85 s; a trivial examination is 2 s wall for 0.7 s of our own time.
Taken apart locally (flat, OneSafe, 315 ms wall): the shell script's per-run
scan of 65 jars for nested libs 82 ms (fixed in ff0208e4: unpacked once into
`plugins-lib/`, atomic rename; the flat run is 230-270 ms since), JVM start to
`main` 120 ms, our run about 95 ms, exit a few ms. Only 2 340 classes load,
some 300 of ours: not classloader hungry. **Twenty of the 65 jars are signed
(Orbit third parties: antlr, aopalliance, javax.activation, ...)** and the JVM
verifies them at every class load, 168 `sun.security` classes: stripping
`META-INF/*.{SF,RSA,DSA,EC}` in a copy of the product took the flat run from
230 to 175 ms (`/data/ythierry/MCC26deploy/flat-unsigned/`). That is the one
cheap packaging win, best done once at install in `install_itstools.sh`
(ITS-Tools-MCC), or after `materialize-products` in the product's pom.
AppCDS (`-XX:ArchiveClassesAtExit`, then `-XX:SharedArchiveFile`) gave 20 ms,
C1 only 15 more: not worth a knob. GraalVM: no Fedora package (Mandrel is a
RHEL build), a CE tarball plus the tracing agent for Guice, Xtext and EMF
would be needed, and a native image runs without C2, so minutes of Java work
slow down to save 0.2 s at start: not for the MCC.
The comparison of the two launchers on a tall node is job 1365883
(`~/MCC26/flat-test/`, `flatbench.sh`, three runs each), queued behind the
campaign; the stderr files were fetched into `warmup-2026-09-07/` for this.
Group quota: files copied to the cluster must take the tree's group, hence
the recipe's `rsync --no-g --chmod=Dg+s`; without it a copy hits "Disk
quota exceeded" at once, the campaign's own logs were never at risk.

## 2026-09-07, 13:55: campaign 2026-09-07 launched (RD, QLA, LTLC, LTLF)

**Product.** The official CI product `202609071104` (ITS-Tools 9fde1d1e, an
empty commit on bae71d95 to rebuild once PetriSpot's `Inv-Linux` carried
3f29014: the first build of bae71d95 had fetched the previous `petri64`;
always check the checksum, `sha256sum .../bin/petri64` against `git show
origin/Inv-Linux:petri64`, here `3843b530652a`). Installed locally with
`install_itstools.sh` in `/data/ythierry/MCC26deploy/MCC-drivers/itstools/`
(the previous local product kept beside it as `itstools.local-bae71d95`),
rsynced with `--delete` over `cluster.lip6.fr:MCC26/MCC-drivers/itstools/`.

**Cluster logs.** The 2026-09-06c results (OS QL QLA RD SM SMA UBA, RDvar)
were copied whole to `/data/ythierry/MCC26archive/2026-09-06c/` (4.7 GB,
file counts checked, README inside) and removed from the cluster, whose
tree is back to its 25 GB of INPUTS.

**Warmup.** All 261 AirplaneLD oracles through `run_oar.sh` (18 instances,
16 examinations), collected in `/data/ythierry/MCC26run/warmup-2026-09-07/`
(`csv/` by the collectors): zero wrong verdict; RC RF RD UB OS QL SM QLA SMA
UBA identical to the campaign tables; **LTLC 288/288 (baseline 287), LTLF
287/288 (baseline 287)** with the wall times down from 240-420 s to 12-22 s
on the COL instances and PT-2000 LTLC from 704 s to 200 s: the glean is back
to nominal. Two Eclipse fatals ("Cannot open display", the shared
`configuration/` race of two jobs of the same instance starting 188 ms
apart): PT-1000 LTLF resubmitted alone answered 16/16, PT-1000 QLA is redone
by the campaign; 2 of 262 is the rate of every campaign so far (30 of 7812
in 06c). CTLC 285/288 and CTLF 273/288 have no baseline of ours (PT-2000
and PT-4000 miss a few within 1800 s).

**Campaign.** `~/MCC26/MCC-drivers/submit-2026-09-07.sh` (a copy is
`Petri/test/mcc/submit-2026-09-07.sh`), detached on the cluster head at
13:52 CEST, log `submit-2026-09-07.log`: RD, QLA, LTLC, LTLF in that order,
one `run_oar.sh` after another, 1800 s / 4 cores / `tall%` /
`runatest_cluster.sh`; the warmup folders of those four moved aside as
`*-before-2026-09-07`. Capacity: 20 tall nodes of 64 cores, about 300
concurrent jobs. Collect each examination when it ends into
`/data/ythierry/MCC26run/2026-09-07/` (rsync without `--delete`, stderr
excluded), run the collectors, point `ITS-Tools latest` at it in
`~/git/MCC-analysis/campaign/example.json`, compare LTLC and LTLF against
`csv/2026-09-06-baseline` (the contest-like ITS-Tools of 2026), then archive
whole and remove from the cluster as above.

## 2026-09-07, 13:30: Effort, the contract of a walker call, and the step cap of a round

**ITS-Tools bae71d95.** `fr.lip6.move.petrispot.runner.Effort { GLEAN, COMMIT }`
is a parameter of every `PetriSpotWalker` call (`runReachability`,
`runBounds`, `runDeadlock`; `runBeside` always commits, the companion's seed)
and of the loop `ReachabilitySolver.applyReductions(reader, doneProps,
timeout, effort)`, which hands it to `randomCheckReachability` and the Parikh
replay. GLEAN caps the sweep at 3 s and the total at 10 s and omits
`--escalate`; COMMIT passes the budgets through with it. Callers: GLEAN for
`AtomicReducer`, `AtomicReducerSR` and `KnowledgeFacts` (the LTL and CTL atoms,
the knowledge loop); COMMIT for the reachability examinations in
`Application`, `GlobalPropertySolver`, `UpperBoundsSolver`, the projection
solver and `DeadlockSolver`. INTEROP.md section 7 and PORTFOLIO's "Three
contracts" record it. The three overloads of `runBounds` and the two of
`runDeadlock` without an effort are gone.

**PetriSpot 3f29014.** The glean alone brought AirplaneLD-PT-0010 LTLC from
209 s to 23 s, but each glean call still spent its whole 10 s: with the
scheduler, a task ending on its `--walkSteps` budget was replaced by a fresh
one until the round's wall time, so a step budget no longer bounded a round
and the driver's "solved nothing on the step budget: stopping" rule never
fired. Now a focused walk carries a step cap of `--walkSteps` per thread over
all its tasks (`Coordinator`, `stepCap`); once spent no task is spawned, the
walk ends with the live ones, `PortfolioResult::stepsExhausted` says so and
`stepBound` reads it. The sweep keeps its clock alone. Effects: Airplane RC
as a glean (`--walkSteps=10000 --sweepTime=3 --totalTime=10`) ends in 4 s
with the stopping message instead of 10 s; as a commit (`--escalate`,
15 s) it now runs 5 rounds raising the budget twice where the old binary ran
one round for the wall time, same 3 witnesses. RERS17pb114 QLA (15 812 vs
16 760 in 4.5 s) and Erlangen full QLA (759 vs 832 in 12 s) are unchanged in
kind: both are settled by the sweep. Through the harness with the new
`petri64` in the deploy tree (`smoke-airplane-{LTLC,RC}-effort2.log`):
**Airplane LTLC 16/16 in 9.5 s** (the two glean calls 3 s each, the 3 s
sweep cap), Airplane RC 16/16 in 50 s as before (one commit call of 50 s on
two properties the walker cannot settle), CTLC 16/16 in 12 s.

Open, seen on Erlangen: with 78 000 open targets after the sweep the driver
finds no round "worth a walk" (under 20 ms each) and hands back 18 of the 30
seconds, `--escalate` or not; PORTFOLIO's G2 (budgets from the clock) covers
it, the sweep should simply go on with the time left.

The deploy tree: `MCC-drivers/itstools/` was wiped by mistake this session
(a `rm -rf` one level too high) and restored by `install_itstools.sh` (a
fresh clone of ITS-Tools-MCC), `install_greatspn.sh` and `mkdir bin`; the
product inside is the local build of bae71d95 (`202609071031`) with
`petri64` replaced by the local 3f29014 build (the CI's copy kept as
`/data/ythierry/MCC26deploy/petri64-1031.bak`). The kept CI product
`itstools.ci-1624` is gone; the cluster still runs that product.

## 2026-09-07, 12:40: inf-stutter, the third spotutil job

`spotutil inf-stutter FILE.hoa` (Spot-BinaryBuilds 6258f61, Linux CI green,
gh-pages binary republished): per state q, the letters x whose infinite word
x x x ... is accepted from q, enumerated (up to 14 atomic propositions, else
status 3) through a one-state word automaton in product with the automaton
restarted at q; the answer is an HOA of self-loops over the same states.
ITS-Tools fcd6e174 `SpotRunner.computeInfStutter` is that one call, the
per-state loop of three Spot processes (`autfilt --small`, `ltl2tgba` of the
stuttering formula, `autfilt --product-and`) and its `simplify` helper are
gone; when the tool fails every state gets `false`, which both consumers
(`RandomProductWalker`, `KnowledgeFacts.computeEGknowledge`) read as nothing
known. Checked on AirplaneLD-PT-0010 LTLC against the previous local log
(`/data/ythierry/MCC26deploy/smoke-airplane-LTLC-{local,infstutter}.log`):
the 15 lists are the same formulas (conjunct order aside), same 16 verdicts,
stuttering time 843 ms to 141 ms; Sudoku-PT-AN01 LTLF 16/16 in 2.8 s. The
deploy tree `itstools/` now holds product `202609071016` with the locally
built (unstripped, 114 MB) spotutil copied over the fetched one; the next
local `mvn` build fetches the published binary, which has the subcommand.
Not yet observed in a run: `stutter-states` and `sensitivity` (no run so far
printed their call). Wall time of Airplane LTLC unchanged at 209 s: the
walker budget (PORTFOLIO "three contracts") remains the pending piece.

## 2026-09-07, 12:15: spotutil wired end to end, the local product in the deploy tree

Spot-BinaryBuilds builds Spot 2.16 and `spotutil` (C++20: the 2.16 headers
use `std::set::contains`; a first CI run compiled it as C++17 and the
gh-pages deploy with `clean` dropped the binary, which made the ITS-Tools
fetch fail). ITS-Tools 1ed32177 (the plugin, `SpotRunner`) and 7c39bc62 (the
`p2.inf` chmod touchpoints, which named the deleted scripts and broke the
product install) are pushed; the CI run of 7c39bc62 was re-run through the
API once the binary was published. Locally: `mvn -o install -DskipTests`
succeeds (`/data/ythierry/itstools-mvn17.log`, 1:26 min), the product
`202609070959` carries `spotutil-linux64` with the three 2.16 Spot binaries
and the master `petri64` (`8de6bbfccc06`, the scheduler and coordinator, not
yet run at scale). **That local product now sits in
`/data/ythierry/MCC26deploy/MCC-drivers/itstools/itstools/`**; the CI product
of 1624 with the 6a5d288 walker is kept beside it as `itstools.ci-1624`.
Through `run_test.pl`: Sudoku-PT-AN01 LTLF (2 tracebacks in the campaign)
16/16 in 4.4 s, no traceback; AirplaneLD-PT-0010 LTLC 16/16, no traceback,
still 210 s on the two walker calls of 135 s and 70 s: the Java budget change
of PORTFOLIO's "three contracts" is the pending piece.

## 2026-09-07: spotutil, the three contracts written

**spotutil** (Spot-BinaryBuilds `tools/`, commits b4a9e6b and after): one
static binary over the Spot C++ API with two subcommands, `stutter-states`
(was `autstates.py`) and `sensitivity` (was `senseclsl.py`), same output
formats; CMake against the Spot install of `build_spot.sh`, CLI11 vendored;
tested locally against the Spot 2.14.5 install of that repo on three
automata. The Spot-BinaryBuilds CI publishes it on gh-pages beside
`ltl2tgba`. That repository now builds the **Spot 2.16 release** (2a74e64)
instead of the 2.14.5.dev snapshot: the knowledge integration ITS-Tools
uses (`--given-formula`, `--given-automaton`, `--given-strategy` with
`auto-small` and `auto-si`, `--product-and`, `--included-in`) is in the
release, its `auto` choice gaining a tie-break on transitions. **ITS-Tools commit 1ed32177 (local, not pushed)** fetches
`bin/spotutil-linux64` in the plugin's pom, drops the two scripts, and
`SpotRunner` calls the subcommands; push it once
`https://github.com/yanntm/Spot-BinaryBuilds/raw/gh-pages/spotutil` answers,
else the ITS-Tools build fails on the fetch. Then the product chain and a
warmup on an LTL instance whose log carried the traceback.

**PORTFOLIO.md** gained "Three contracts for a walker call": glean, commit,
companion, with the Java call sites and their contract; no code yet.

## 2026-09-07, morning: LTL slowed by the walker's budget, the Spot scripts, the total pages

**LTL.** `LTLC/OAR.1337407` (AirplaneLD-PT-0010, 15 s in the contest) took
215 s: two walker calls of 135 s and 70 s. They come from
`ReachabilitySolver.randomCheckReachability`, whose budget is
`30 + 5 min(|tocheck|, 50)` seconds and whose `runReachability` call always
passes `--escalate` (added for the deadlock and reachability examinations,
where a walk that returned after one second was the loss): under LTL the
reductions (`applyReductions` with `ReductionType.LTL`/`SI_LTL`,
`Application.java` 717) call it to glean easy atoms, and since 202609060003
the walk spends the whole budget whenever one atom stays open. The
`AtomicReducer`s themselves ask for 30 s. Fix, on the Java side (not made
from here): `randomCheckReachability` takes the `ReductionType`; for LTL,
SI_LTL and CTL it asks the walker for a short glean (`sweepSeconds` 1, total
a few seconds) and no escalation, `PetriSpotWalker.runReachability` gaining
an `escalate` flag; the reachability examinations keep the formula.

**Spot.** 1 108 LTLF logs carry a python traceback: `autstates.py` and
`senseclsl.py` in `fr.lip6.ltl.spot.binaries` do `import spot`, absent on the
cluster's python. They compute stutter-invariant states made forward closed
(`spot.stutter_invariant_states`, `make_stutter_invariant_forward_closed_inplace`,
HOA out) and stutter, lengthening and shortening insensitivity (`closure`,
`sl`, `complement`, `product`, `is_empty`). Both are a page of the Spot C++
API each; built static like `ltl2tgba-linux64` in the same plugin, they would
end the python dependency. Preferred over shipping the python module.

**Total pages** (MCC-analysis, pushed): the completion-against-time plot is
gone (every run at the wall sat on one line), replaced by A against B scatters
of wall time and of atoms answered on the shared instances, an
answered-against-atoms plot per set, the progression along a family on wall
time by default (on completion two complete runs coincide, which was the
"identical curves"), and A/B filters of the runs table.

## 2026-09-07, 03:00: the rerun is complete and collected

All seven examinations of `submit-2026-09-06c.sh` are in
`/data/ythierry/MCC26run/2026-09-06c/`, the tables in
`Petri/test/mcc/csv/2026-09-06c/` with a README (the reading), the pages
rebuilt (`ITS-Tools latest` is that set). The LTLF baseline is complete
(1 953), its tables with RC, RF, LTLC in `csv/2026-09-06-baseline/`, the
whole logs archived to `MCC26archive/2026-09-06-baseline/` and LTLF flushed
from the cluster, whose tree is now INPUTS plus the seven rerun directories.
The cluster queue is empty; no watch is armed. Next campaign: rebuild the
product through the chain for master (the scheduler, the coordinator, the
quest tool: 7acb697 to 243adaa), warm up, then the totals first.

## 2026-09-07, 01:30: Walker.h split, the quest tool, seats

ca706a8 moves the target index to `walk/TargetIndex.h` (Walker.h 474
lines). The next commit is the first cut of step 3b: a `quest` tool (one
target per task, `QuestSweep::pick` at spawn, the task ends on the claim),
spawns built outside the scheduler's lock, a spawn attempt whenever a task
ends, a seat per kind; the auto sweep runs sync, quest and rare under the
shares. Yardsticks, 4 runners: Stigmergy 2 357 (1 554 before), ResIsolation
1 000 in 12 s (8.9 s), Erlangen 1 866 (2 458), RERS 29 042, DLCflexbar
76 160; the trade-off is written in WALK_PLAN 10.12 step 3b. Peterson (the
user asked): QLA answered in full up to PT-4 (690 atoms in 1 216 s), PT-5
954 of 1 242, PT-6 1 698 of 2 030, PT-7 2 535 of 3 096 in the rerun (the
morning: 928, 1 297, 1 887); SMA in full up to PT-4, then 513 of 834, 808 of
1 330, 825 of 1 992 (morning 596, 657, 761): the walker carries most of it
(`walk solved` 1 655 on PT-5 QLA against 954 answered, the rest duplicates of
the identical-property reduction), the diagrams the small instances.

## 2026-09-07, 00:30: the coordinator (step 3a), SMA collected

Commits 42ee6df (WALK_PLAN 10.12, the coordinator design) and ee14593
(`walk/Coordinator.h`): tasks that stop progressing are parked (best state
to the pool, yield to the arm table, walker freed), new tasks are spawned
from (pooled state, tool) pairs by decayed yield, untried pairs first,
nothing spawned once nothing is left. Yardsticks, 4 runners: Erlangen full
QLA 2 417 (all-quest 2 792; equal shares 1 444), Stigmergy 1 609 (best),
RERS 30 186 in 10 s (all-rarity 30 121), ResIsolation 1 000 in 9.6 s,
DLCflexbar 76 160 of 76 160, controls unchanged; CAN gathering still open
after 263 tasks over 139 pairs (the case for the LP refinements). Two policy
bugs found and fixed on the way: fresh tasks rewarded for their own first
firings (novelty is now what nobody fired before) and spawning after every
target was claimed (the `more` hook). Options `--tasks` (live cap), `--grant`,
`--spawn`, `--shares`, `--shareFloor`, `--slice`, `--sliceMs`.

Housekeeping owed: `Walker.h` is 541 lines; the target index (own, up and
down lists, the checks) should move to `TargetIndex.h`. Step 3b (quests as
spawn decisions, child tasks with a budget) and step 4 (LP tasks) follow.

**SMA rerun collected** into `2026-09-06c/SMA` (1 681 finished, tables in
`csv/2026-09-06c/`, pages rebuilt): 59 models better, 42 worse, answered
4.25 M -> 3.81 M; most of these runs used the 1624 binary (the swap to
6a5d288 was at 21:25, SMA ran 19:20 to 22:17), so read them with the QLA
caveats. UBA, RD, QL, SM, OS still running.

## 2026-09-06, 23:00: time sharing, step 1 landed

WALK_PLAN.md 10.11 (rewritten as agreed: a pool of exploration tasks over
runner threads, shares that follow results, subquests as child tasks with a
budget) and PORTFOLIO.md's new section are the design; commit 7acb697 is
step 1: `walk/Task.h`, `WalkTask.h`, `Scheduler.h`, `Walker::begin/runSlice/
finish`, options `--tasks --slice --sliceMs`, a report line per task and per
strategy kind (steps, running ms, steps/ms, claims, claims/s, slices, capped
slices). Equal shares, so with twice the tasks each kind gets half the time.
Parity with the threads held on every control; the yardsticks and the
calibration figures are in 10.11 step 1 and `Petri/test/logs/sched-*.log`.
Two figures to keep in mind for step 2 (adaptive shares) and the sampling
question: on Erlangen the rarity tasks pay 50 000 arc visits and 25 000
target checks a step (0.4 steps/ms), on RERS the quest tasks fire a cheap
local loop at 130 steps/ms with 1 arc visit a step and claim nothing.
Walker.h is 521 lines; the target index is the piece to extract next.

**Step 2 landed too** (44a1eb0): shares follow the claims per running
second, decayed over three seconds, above a floor of a tenth; two tasks per
runner by default. One configuration now reaches RERS 29 352 (rarity 91 %),
ResIsolation 1 000 targets in 10.6 s (quests 86 %), Erlangen full QLA 2 140
(quests 94 %), DLCflexbar 76 160 of 76 160. Steps 3 (subquests as child tasks
with a budget) and 4 (LP tasks) are next; the binary on the cluster is
6a5d288, before the scheduler.

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

**The rerun in progress (SMA, UBA, RD, QL, SM, OS) started on product
`202609061624`, which carries both problems.** At 21:25 the `petri64` inside
the deployed product was replaced, without an ITS-Tools rebuild, by the
CI-built static binary of master `6a5d288` taken from `Inv-Linux` (checksum
`c90125eb9772`; the 1624 binary kept as
`/data/ythierry/MCC26deploy/product-check/petri64-1624`), and the tree rsynced
to the cluster: jobs starting after that, and running jobs at their next
walker call, use the fixed walker; QLA and the SMA/UBA runs before 21:25 do
not. A log's binary is told by its `Version` line only for the Java side;
the walker's thread report line says `condemned` from `5ade6af` on and the
`Components: N of M semiflows` line from `6a5d288` on. Local models for the cases: `/data/ythierry/rers/`
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
