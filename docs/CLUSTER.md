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
| `/data/ythierry/MCC26deploy/MCC-drivers/` | the local deploy tree, built here and rsynced up, one subtree at a time |
| `/data/ythierry/MCC26run/<date>/<EXAM>/` | collected logs, `csv/` beside them |
| `/data/ythierry/MCC26archive/<date>/` | finished campaigns, moved out of `MCC26run` |

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

The image is compiled for AVX2, so it runs on `tall%` only. On `small%` it
prints `The current machine does not support all of the following CPU features`
and exits having answered nothing -- a silent zero-`FORMULA` log, not a crash.
Until `build-native.sh` gets an `-march`, a campaign on the image is `tall%`.

What actually ran is in the job's `.stderr` (`runeclipse.sh` traces with
`set -x`), not in the log:

```
grep -m1 'exec .*its-tools-native' <EXAM>/OAR.<jobid>.stderr    # native, else the Eclipse launcher
```

Startup on `tall11`, OneSafe on AirplaneLD-PT-0010, cold then warm
(`~/MCC26/flat-test/flatbench.sh`): Eclipse 3094 then ~1000 ms, flat 930 then
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
is the old nodes, no AVX2). `TAG=<name>` qualifies the result folder
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
down). For several examinations, one `run_oar.sh` after another in a detached
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

## 4. Collect, as often as wanted

```
BASELINE=Petri/test/mcc/csv/<run to compare with> bash Petri/test/mcc/collect.sh <date> QLA LTLC LTLF
```

Rsyncs the folders into `/data/ythierry/MCC26run/<date>/` (never with
`--delete`; a draining campaign is collected as it goes and a rerun adds the
new logs), runs `mcclogs2csv.py` on the classic examinations and
`totallogs2csv.py` on `QLA SMA UBA`, then `report.py` into
`<date>/csv/REPORT.md`, the readable result, against `$BASELINE` when set.
`toolsupport.py` (which contest tools back a value we miss) is the manual
extra:

```
python3 Petri/test/mcc/toolsupport.py ~/git/pnmcc-models-2026/website/raw-result-analysis.csv /data/ythierry/MCC26run/<date>/csv/verdicts.csv -o /data/ythierry/MCC26run/<date>/csv/support.csv
```

The tables of a campaign worth keeping go to `Petri/test/mcc/csv/<date>/`
with a README saying what it showed. The Eclipse fatals (about 5 per
examination) sit in the `failure` column: read it before believing a one
second miss. Run to run variance is a few instances per examination; no A/B
worth less than a few dozen instances can be read from one run.

## 5. Pages, then archive

The pages (`~/git/MCC-analysis/campaign/`, its README) are built from
`campaign/example.json`: one *set* per campaign, `"logs":
["/data/ythierry/MCC26run/<date>/*"]`, so a rebuild picks up whatever the
folder holds. A new campaign date is a new set at the top of the list (copy
the previous block, rename it); the previous `latest` keeps its dated name.

```
bash Petri/test/mcc/collect.sh <date> QLA LTLC LTLF --pages     # collect and rebuild the pages (about four minutes)
python3 ~/git/MCC-analysis/campaign/build.py ~/git/MCC-analysis/campaign/example.json   # the rebuild alone
python3 ~/git/MCC-analysis/campaign/serve.py /data/ythierry/MCC26run/pages --port 8080  # browse, logs served from disk
```

The order is: collect, update `example.json` if the campaign is new, rebuild
the pages, browse. When the campaign is complete and read, archive it and free
the cluster. The cluster home is under quota and the inputs alone are 25 GB
that do not compress; result folders are deleted there as soon as they are
archived, `/data` has terabytes and keeps everything.

```
mv /data/ythierry/MCC26run/<date> /data/ythierry/MCC26archive/<date>      # then point the set's glob at the archive
ssh cluster.lip6.fr 'cd MCC26/MCC-drivers && rm -rf QLA LTLC LTLF'        # only after the archive is complete
ssh cluster.lip6.fr 'du -sh MCC26; quota -s 2>/dev/null | tail -2'         # what is left up there
```

## 6. Run one instance locally through the same harness

```
cd /data/ythierry/MCC26deploy/MCC-drivers
BK_TOOL=itstools ./run_test.pl oracle/AirplaneLD-PT-0010-LTLC.out -t 300      # one examination
BK_TOOL=itstools ./run_model.sh AirplaneLD-PT-0010 -t 60                       # every examination with an oracle
BK_TOOL=petrispot PETRISPOT_BIN_RUN=~/git/PetriSpot/build/petri64 ./run_test.pl oracle/AirplaneLD-PT-0010-RC.out -t 60
```

The log is the same as a cluster log, so the collectors read it too.
