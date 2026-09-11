# Guide: the cluster, day to day

How a campaign is run, watched, collected and archived. The harness itself
(what the tree holds, how it is built and rsynced, the OAR details) is in
`BENCH.md`; this is the operating sequence. Every command is meant to be run
as written.

The head node `cluster.lip6.fr` is for `oarsub` and for deployment only. It
has no tools beyond the basics (no Java, no compiler, `zip` and `rsync` about
sum it up), and it is not meant for test runs, not even a small one: anything
that executes a tool is a job. Deployment works because the compute nodes see
the same file system as the head we log into: what is rsynced to
`~/MCC26/MCC-drivers/` is what every job finds, and the result folders the
jobs write are read back from the head.

## Layout

| where | what |
| --- | --- |
| `cluster.lip6.fr:~/MCC26/MCC-drivers/` | the harness: `run_oar.sh`, `oracle/`, `INPUTS/`, one folder per tool (`itstools/`, `petrispot/`, `hsc/`...) |
| `.../MCC-drivers/<EXAM>/` | the results of one examination: `OAR.<jobid>.stdout` (the log) and `.stderr` (the `time -p` trailer); they exist only on the cluster until collected. `<EXAM>.<tag>/` when the submission carried `TAG` (one tool, or one setting, per folder) |
| `cluster.lip6.fr:~/MCC26/flat-test/` | the launcher bench (Eclipse, flat, native image) |
| `/data/ythierry/MCC26deploy/MCC-drivers/` | the local deploy tree, built here and rsynced up, one subtree at a time; its `oracle/` is the only oracle copy outside `pnmcc-models-2026` |
| `/data/ythierry/MCC26deploy/native-v2/its-tools-native-v2` | the x86-64-v2 image for `small%`; `hsc-sweep/` the libHSC order sweep's deploy; `pnmcc-tests-src/` the test harness sources |
| `/data/ythierry/MCC26logs/<tool>/<build>/<EXAM>/` | collected logs, `csv/` beside them (its `README.md` is the index); `web/campaign/` and `web/order-sweep/` the generated pages, `local/<name>/` a local experiment's files, `total/` the total-examination tables |
| `/data/ythierry/scratch/` | one-shot tests, emptied at will |

Nothing is written at the top of `/data/ythierry` (its `README.md`,
`CLAUDE.md`): no loose log, no unpacked model, no oracle copy. A model is
unpacked by the harness under `INPUTS/` or by hand under the repo's
git-ignored `bench/models/`, and deleted after use anywhere else.

Examination folder names: `RC RF RD UB L QL SM OS SS LTLC LTLF CTLC CTLF`
and the total examinations `QLA SMA UBA` (`TOTAL_QUERIES.md`).

## 1. Deploy a tool

The cluster never builds. Install locally, then rsync the subtree (never the
root once results exist: `--delete` would erase them).

```
cd /data/ythierry/MCC26deploy/MCC-drivers
rm -rf itstools && ./install_itstools.sh && (cd itstools && ./install.sh)   # ITS-Tools: the published product + its-tools-native
(cd petrispot && ./install.sh)                                              # petri64 from the PetriSpot CI; PETRISPOT_BIN=<path> for a local build
rsync -rlptD --no-g --chmod=Dg+s --delete /data/ythierry/MCC26deploy/MCC-drivers/itstools/ cluster.lip6.fr:MCC26/MCC-drivers/itstools/
```

`--chmod=Dg+s` and `--no-g` are load bearing (the quota group, `BENCH.md`).
Do not touch a tool folder while a campaign that uses it is running.

### Which launcher a campaign gets

`install_itstools.sh` downloads the product **and** `its-tools-native`, and
`runeclipse.sh` execs the native image whenever that file is present. A fresh
install therefore switches the launcher of the next campaign without saying so.
Decide it deliberately: keep the file for the native image, delete it from the
deploy tree before the rsync for the Eclipse launcher.

The CI image is compiled for x86-64-v3 (native-image's default, AVX2), so it
runs on `tall%` only. On `small%` and `big%` it prints `The current machine
does not support all of the following CPU features` and exits having answered
nothing -- a silent zero-`FORMULA` log, not a crash. An image for those nodes
is built locally, `NATIVE_MARCH=x86-64-v2` to `build-native.sh` on the deployed
product (40 s, 10 GB), and put in the deploy tree as `its-tools-native`
before the rsync; it runs everywhere (v2 is a subset of v3). Tested on
`small10` and `big12` on AirplaneLD-PT-0010, OneSafe, CTLCardinality and
ReachabilityCardinality answering, the external binaries found through
`plugins/`. Nothing else in the tree minds the nodes: `petri64`, `hsc-pn` (the
CI's static GMP is Ubuntu's fat build) and the Eclipse launcher on the nodes'
Java 21 all run on `small%`.

### The node classes, and what a job gets

Probed with `Petri/test/mcc/probe_node.sh` (one job per class and core count):

| class | nodes | CPU | cores (threads) | RAM | x86-64 | RAM per core |
| --- | --- | --- | --- | --- | --- | --- |
| `small%` | 24 | 2 x Xeon E5645, 2.4 GHz (Westmere, 2010) | 24 (12 physical, HTT) | 62 GiB | v2, no AVX | 2.6 GiB |
| `big%` | 26 | 2 x Xeon X5690, 3.47 GHz (Westmere) | 24 (12 physical, HTT) | 141 GiB | v2, no AVX | 5.9 GiB |
| `tall%` | 20 | current, AVX2 | 64 | | v3 | |

OAR grants no CPU time limit but **a memory cap proportional to the cores
requested**: `memory.max` on the job's systemd slice
(`/sys/fs/cgroup/oar.slice/oar-<uid>.slice/oar-<uid>-j<job>.slice`) is the
node's RAM divided by its cores, times the cores of the job; the shell's own
scope shows `max`, the slice above it holds the number. Measured on `small10`:
3 cores 7.9 GiB, 6 cores 15.7 GiB; 1 core cannot hold 4 GB, 12 cores hold 28.
So a job that needs the 16 GB of an MCC run asks for **6 cores on `small%`**
(a quarter of the node, 3 physical cores) and 3 on `big%`; a core is a
hyper-thread there, so the count is not the parallelism it reads as.
The `small` nodes sit in standby (`oarnodes`: `Absent (standby)`) and wake for
a job in about a minute.

What actually ran is in the job's `.stderr` (`runeclipse.sh` traces with
`set -x`), not in the log:

```
grep -m1 'exec .*its-tools-native' <EXAM>/OAR.<jobid>.stderr    # native, else the Eclipse launcher
```

Startup on `tall11`, OneSafe on AirplaneLD-PT-0010, cold then warm
(`Petri/test/mcc/flatbench.sh`, run on a node): Eclipse 3094 then ~1000 ms, flat 930 then
~620 ms, native 596 then ~25 ms; a bare JVM starts in 45 ms. Against an 1800 s
timeout the saving is noise; the image is worth it for short examinations.

The image carries a closed world: a class reached by reflection that the
tracing agent never recorded is a `MissingReflectionRegistrationError` at run
time, and the engine that hit it answers nothing while the rest of the run
looks healthy. The recipe (trace the failing run, rebuild) is in
`~/git/ITStools/ITS-commandline/native/README.md`.

## 2. Warm up, then submit

Always one small instance first: it catches a broken install in a minute.

```
ssh cluster.lip6.fr
cd ~/MCC26/MCC-drivers
BK_TOOL=itstools ./run_oar.sh 'oracle/AirplaneLD-PT-0010-*.out'                    # warmup, all examinations of one instance
TIMEOUT=1800 WALLTIME=0:45:0 CORES=4 HOSTS="tall%" BK_TOOL=itstools ./run_oar.sh 'oracle/*-RD.out'   # one examination, 1953 jobs
```

`TIMEOUT` is the budget `run_test.pl` gives the tool, `WALLTIME` the OAR
limit (comfortably above), `HOSTS` `tall%` for the current hardware (`small%`
and `big%` are the old nodes, x86-64-v2: the node table above says what runs
there and how many cores 16 GB take). `TAG=<name>` qualifies the result folder
(`SS.hsc`, `RC.itstools`), so two tools, or two settings of one tool, run
the same examination side by side without mixing logs; `collect.sh` takes
the qualified name as its examination argument and the collectors read the
examination from the log, not from the folder:

```
TIMEOUT=300 WALLTIME=0:10:0 CORES=4 HOSTS="tall%" TAG=itstools BK_TOOL=itstools ./run_oar.sh 'oracle/*-SS.out'
bash Petri/test/mcc/collect.sh <date> SS.itstools   # the folder name, tag included
```
 Submit from the foreground and let it end (a few
hundred `oarsub` take minutes; parallel submissions have brought the head
down). **Never rewrite `run_oar.sh` (or any script a job or a loop is
reading) while it runs**: bash reads a script as it goes, so an overwrite
mid-loop makes it die on a syntax error at the changed offset, silently
truncating a submission. Copy a new version up once the loop has ended, and
check what actually got submitted (the `comm -3` recipe of `BENCH.md`,
oracle list against the jobs in `oarstat -f -u` plus the logs already
written). For several examinations, one `run_oar.sh` after another in a detached
script, as `~/MCC26/MCC-drivers/submit-<date>.sh` does (a copy of the last
one is `Petri/test/mcc/submit-2026-09-07.sh`).

## 3. Watch

```
bash Petri/test/mcc/cluster_status.sh              # every result folder; or name examinations: cluster_status.sh QLA LTLC
```

One line per examination folder: logs written (one `OAR.<id>.stdout` per
finished job) and how many carry the harness trailer (a job killed at the
wall leaves none), then the jobs by state (`R` running, `W` waiting). The
18-log folders are the Airplane warmups. About 450 jobs an hour at 4 cores on
`tall%`. Resubmit only the killed ones (`BENCH.md`, the `comm -3` recipe).

## 4. Reap: collect into the one log root

Every log lives under **`/data/ythierry/MCC26logs/<tool>/<build>/<EXAM>[.tag]/`**
and nowhere else (its `README.md` is the rule and the index). `<tool>` is the
harness's tool folder (`itstools`, `hsc`, `petrispot`); `<build>` is the
**product build id**, the campaign's label — for ITS-Tools the plugin
timestamp of what is deployed, read on the head node before submitting:

```
ssh cluster.lip6.fr 'ls MCC26/MCC-drivers/itstools/plugins | grep -o "1.0.0.20[0-9]*" | head -1'
```

(a tool without a build id in its logs is labelled by date, `hsc/20260908`).
A campaign that mixed builds keeps the one grabbed: old logs are not rerun,
they stay informative.

```
BASELINE=/data/ythierry/MCC26logs/itstools/<previous build>/csv bash Petri/test/mcc/collect.sh itstools/<build> QLA LTLC LTLF
```

Rsyncs the folders into `MCC26logs/itstools/<build>/` (never with `--delete`;
a draining campaign is collected as it goes and a rerun adds the new logs),
runs `mcclogs2csv.py` on the classic examinations and `totallogs2csv.py` on
`QLA SMA UBA`, then `report.py` into `<build>/csv/REPORT.md`, the readable
result, against `$BASELINE` when set. The same build running the same
examination twice: rename the first folder with a tag before collecting
again (`mv RD RD.precampaign`). `toolsupport.py` (which contest tools back a
value we miss) is the manual extra:

```
python3 Petri/test/mcc/toolsupport.py ~/git/pnmcc-models-2026/website/raw-result-analysis.csv /data/ythierry/MCC26logs/itstools/<build>/csv/verdicts.csv -o /data/ythierry/MCC26logs/itstools/<build>/csv/support.csv
```

The tables of a campaign worth keeping go to `Petri/test/mcc/csv/<build>/`
with a README saying what it showed. The Eclipse fatals (about 5 per
examination) sit in the `failure` column: read it before believing a one
second miss. Run to run variance is a few instances per examination; no A/B
worth less than a few dozen instances can be read from one run.

## 5. Pages, README, then free the cluster

The pages (`~/git/MCC-analysis/campaign/`, its README) are built from
`campaign/example.json`: one *set* per campaign, `"logs":
["/data/ythierry/MCC26logs/itstools/<build>/*"]`; the oracle is the deploy
tree's (`MCC26deploy/MCC-drivers/oracle`), the pages are built into
`MCC26logs/web/campaign/`. A new campaign is a new set at the
top of the list (copy the previous block, rename it).

```
bash Petri/test/mcc/collect.sh itstools/<build> QLA LTLC LTLF --pages   # collect and rebuild the pages (about four minutes)
python3 ~/git/MCC-analysis/campaign/build.py ~/git/MCC-analysis/campaign/example.json   # the rebuild alone
python3 ~/git/MCC-analysis/campaign/serve.py            # browse http://127.0.0.1:8080/campaign/ and /order-sweep/, logs served from disk
```

When the campaign is complete and read: write `MCC26logs/<tool>/<build>/README.md`
(what ran, the timeout, what it showed, the wrong verdicts), add its line to
the table in `MCC26logs/README.md`, then free the cluster — the cluster home
is under quota and the inputs alone are 25 GB that do not compress; result
folders are deleted there as soon as they are collected in full (compare the
log counts first), `/data` has terabytes and keeps everything.

```
ssh cluster.lip6.fr 'cd MCC26/MCC-drivers && rm -rf QLA LTLC LTLF'        # only after the counts match
ssh cluster.lip6.fr 'du -sh MCC26; quota -s 2>/dev/null | tail -2'         # what is left up there
```

The former `MCC26run` / `MCC26archive` split was folded into this root on
2026-09-09 (`Petri/test/mcc/reap_reorg.sh`).

## 6. Run one instance locally through the same harness

```
cd /data/ythierry/MCC26deploy/MCC-drivers
BK_TOOL=itstools ./run_test.pl oracle/AirplaneLD-PT-0010-LTLC.out -t 300      # one examination
BK_TOOL=itstools ./run_model.sh AirplaneLD-PT-0010 -t 60                       # every examination with an oracle
BK_TOOL=petrispot PETRISPOT_BIN_RUN=~/git/PetriSpot/build/petri64 ./run_test.pl oracle/AirplaneLD-PT-0010-RC.out -t 60
```

The log is the same as a cluster log, so the collectors read it too.
