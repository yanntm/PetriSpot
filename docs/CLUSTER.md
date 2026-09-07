# Guide: the cluster, day to day

How a campaign is run, watched, collected and archived. The harness itself
(what the tree holds, how it is built and rsynced, the OAR details) is in
`BENCH.md`; this is the operating sequence. Every command is meant to be run
as written; the head node `cluster.lip6.fr` runs nothing but `oarsub`.

## Layout

| where | what |
| --- | --- |
| `cluster.lip6.fr:~/MCC26/MCC-drivers/` | the harness: `run_oar.sh`, `oracle/`, `INPUTS/`, one folder per tool (`itstools/`, `petrispot/`, `hsc/`...) |
| `.../MCC-drivers/<EXAM>/` | the results of one examination: `OAR.<jobid>.stdout` (the log) and `.stderr` (the `time -p` trailer); they exist only on the cluster until collected |
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
is the old nodes, no AVX2). Submit from the foreground and let it end (a few
hundred `oarsub` take minutes; parallel submissions have brought the head
down). For several examinations, one `run_oar.sh` after another in a detached
script, as `~/MCC26/MCC-drivers/submit-<date>.sh` does (a copy of the last
one is `Petri/test/mcc/submit-2026-09-07.sh`).

## 3. Watch

```
ssh cluster.lip6.fr 'cd MCC26/MCC-drivers && for d in RD QLA LTLC LTLF; do echo "$d: $(ls $d 2>/dev/null | grep -c stdout) logs"; done; oarstat -u | tail -n +3 | awk "{print \$2}" | sort | uniq -c'
```

One `.stdout` per finished job; `R` running, `W` waiting. About 450 jobs an
hour at 4 cores on `tall%`. A job killed at the wall leaves a log without a
`Test :` trailer; resubmit only those (`BENCH.md`, the `comm -3` recipe).

## 4. Collect

```
rsync -rz --exclude='*.stderr' cluster.lip6.fr:MCC26/MCC-drivers/{RD,QLA,LTLC,LTLF} /data/ythierry/MCC26run/<date>/
cd /data/ythierry/MCC26run/<date>
python3 ~/git/PetriSpot/Petri/test/mcc/mcclogs2csv.py RD LTLC LTLF -o csv                 # classic examinations
python3 ~/git/PetriSpot/Petri/test/mcc/totallogs2csv.py QLA SMA UBA -o csv --oracles csv/oracles   # total examinations
python3 ~/git/PetriSpot/Petri/test/mcc/report.py csv --baseline <csv folder of the run to compare with>
python3 ~/git/PetriSpot/Petri/test/mcc/toolsupport.py ~/git/pnmcc-models-2026/website/raw-result-analysis.csv csv/verdicts.csv -o csv/support.csv
```

`REPORT.md` in `csv/` is the readable result; the tables of a campaign worth
keeping go to `Petri/test/mcc/csv/<date>/` with a README saying what it
showed. The Eclipse fatals (about 5 per examination) sit in the `failure`
column: read it before believing a one second miss. Run to run variance is a
few instances per examination; no A/B worth less than a few dozen instances
can be read from one run.

## 5. Publish and archive

```
# the pages: point "ITS-Tools latest" at the new folder, name the previous one after its date
$EDITOR ~/git/MCC-analysis/campaign/example.json && (cd ~/git/MCC-analysis && python3 build.py ...)   # see its README
mv /data/ythierry/MCC26run/<date> /data/ythierry/MCC26archive/<date>
ssh cluster.lip6.fr 'cd MCC26/MCC-drivers && rm -rf RD QLA LTLC LTLF'        # only after the archive is complete
```

## 6. Run one instance locally through the same harness

```
cd /data/ythierry/MCC26deploy/MCC-drivers
BK_TOOL=itstools ./run_test.pl oracle/AirplaneLD-PT-0010-LTLC.out -t 300      # one examination
BK_TOOL=itstools ./run_model.sh AirplaneLD-PT-0010 -t 60                       # every examination with an oracle
BK_TOOL=petrispot PETRISPOT_BIN_RUN=~/git/PetriSpot/build/petri64 ./run_test.pl oracle/AirplaneLD-PT-0010-RC.out -t 60
```

The log is the same as a cluster log, so the collectors read it too.
