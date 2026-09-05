# Benchmarking PetriSpot on the MCC 2026 models

PetriSpot is benchmarked as it is used: through ITS-Tools, which carries the
`petri64` binary in its `fr.lip6.petrispot.binaries` plugin and calls it from
the reachability, upper-bounds and structural-reduction paths. The measurement
harness is the MCC one, so the numbers are directly comparable with the
contest: `run_test.pl` of `pnmcc-tests` drives `BenchKit_head.sh` of
`MCC-drivers` on one model instance and one examination, and compares every
verdict against the 2026 oracle.

## The four repositories

| repository | role |
| --- | --- |
| `MCC-drivers` | the `BenchKit_head.sh` dispatcher and one folder per tool; `BK_TOOL` selects the folder, the `xred` suffix runs the ITS-Tools reducer first |
| `pnmcc-tests` | `run_test.pl`, `runatest.sh`, `runatest_cluster.sh`, `limit_time.pl`, `analysis/` |
| `pnmcc-models-2026` | the model archives (`website/INPUTS/*.tgz`) and the consensus oracles (`website/oracle.tar.gz`), published on its `gh-pages` branch |
| ITS-Tools | not cloned: the linux product is taken from the CI at `https://lip6.github.io/ITSTools/` |

Everything binary comes from a CI artifact, never from a local build: the
ITS-Tools product from the ITSTools CI, `petrispot/bin/petri64` from the
PetriSpot CI (`Inv-Linux` branch), the other tools from their own installers.
So a benchmark measures a published state of the tools, reproducible from the
same URLs. To benchmark a PetriSpot change, publish it on the PetriSpot CI and
let the ITS-Tools CI pick it up.

## Building the tree

Built locally in `/data/ythierry/MCC26deploy/`, then rsynced to
`cluster.lip6.fr:~/MCC26/`. The cluster head runs nothing: it only holds the
files and submits `oarsub` jobs.

```
DEPLOY=/data/ythierry/MCC26deploy
mkdir -p $DEPLOY && cd $DEPLOY
git clone https://github.com/yanntm/MCC-drivers.git MCC-drivers
git clone https://github.com/yanntm/pnmcc-tests.git pnmcc-tests-src

cd MCC-drivers
# the harness must sit next to BenchKit_head.sh (README of pnmcc-tests)
cp ../pnmcc-tests-src/{run_test.pl,runatest.sh,runatest_cluster.sh,limit_time.pl,install_input.sh,install_oracle.sh} .
cp -r ../pnmcc-tests-src/analysis .

./install_itstools.sh    # clones ITS-Tools-MCC into itstools/
for i in */ ; do [ -f $i/install.sh ] && (cd $i && ./install.sh) ; done

# models: pre-seeded so no node ever needs the network. install_input.sh only
# downloads when INPUTS/<model>.tgz is absent.
mkdir -p INPUTS && cp ~/git/pnmcc-models-2026/website/INPUTS/*.tgz INPUTS/
```

The oracles are committed in `MCC-drivers/oracle/` (24827 files, `ORACLE2026`
and `TEDD2026` verdicts), so `install_oracle.sh` is not needed. The oracles of
the total examinations (`TOTAL_QUERIES.md`) are not committed, they are 400 MB
of `?`: `./install_total_oracles.sh` fetches them into the same folder.

`INPUTS/` holds the whole 2026 set. Two archives, `StigmergyCommit-PT-11b`
(140 MB) and `TokenRing-PT-050` (99 MB), exceed the GitHub pages file size
limit, so `install_inputs.sh` deletes them and `pnmcc-models-2026` does not
publish them: a tool that is handed the model folder finds no `model.pnml` and
the whole test reports nothing. They come from the contest archive instead,
which carries every model at full size:

```
cd ~/git/pnmcc-models-2026/website
wget https://mcc.lip6.fr/2026/archives/INPUTS-2026.tar.gz          # 1.1 GB
tar xzf INPUTS-2026.tar.gz INPUTS-2026/StigmergyCommit-PT-11b.tgz \
                           INPUTS-2026/TokenRing-PT-050.tgz
```

Unpack, run `patch_models.pl` on the folder and repack, as `install_inputs.sh`
does for every other model, then copy the two `.tgz` into `INPUTS/`. The models
themselves are byte identical from one edition to the next; the property files
are not, so the archive must be the 2026 one.

## Rsync to the cluster

```
rsync -rlptD --no-g --chmod=Dg+s --delete /data/ythierry/MCC26deploy/ cluster.lip6.fr:MCC26/
```

`--delete` over the whole tree is destructive once a campaign has run: the
result directories `L/`, `RD/`, `UB/` and their `OAR.*.stdout` exist only on
the cluster. Rsync a subtree (`.../itstools/`, `.../INPUTS/`) rather than the
root when results are in place.

Until ITS-Tools `202609052009` the product zip stored every plugin binary as
mode `644` -- `petri64`, `its-reach-linux64`, `louvain-linux64` alike -- so
after an install only the Java runtime ever made them executable, and a hundred
concurrent jobs sharing one tree raced on that. A `p2.inf` chmod touchpoint in
each binaries plugin now sets 755 when p2 materializes the product, and a fresh
install comes out executable. On an older product, do it by hand before the
rsync:

```
chmod +x MCC-drivers/itstools/itstools/plugins/*/bin/*
```

A job that cannot spawn `petri64` says so once
(`PetriSpot I/O error (PFLOWS): ... Exec failed`) and then carries on without
PetriSpot, so the loss is silent in the verdicts.

`--chmod=Dg+s` is load bearing, and `-a` is wrong here. The home directory on
the cluster is setgid to the group `ythierry`, which is the group that carries
the disk quota; the user's primary group `MoVe` has none. A plain `rsync -a`
copies the local directory modes over, dropping the setgid bit, after which
every file created under the tree is charged to `MoVe` and the write fails with
`Disk quota exceeded` -- including `mkdir`, while `touch` of an empty file still
works, since an empty file costs no block. If that happens, the repair is
`chgrp ythierry <dir> && chmod g+s <dir>`.

## Running

`run_model.sh` (in the tree, next to `run_test.pl`) runs every examination we
have an oracle for on one instance:

```
BK_TOOL=itstools ./run_model.sh AirplaneLD-PT-0010 -t 60
```

A single examination is the plain `pnmcc-tests` call:

```
BK_TOOL=itstools ./run_test.pl oracle/AirplaneLD-PT-0010-RF.out -t 300
```

On the cluster, use `runatest_cluster.sh` in place of `runatest.sh`: it
unpacks the model into a `$$`-suffixed folder so concurrent jobs on the same
model do not share a working directory.

`run_oar.sh` submits the campaign: one `oarsub` per oracle file, that is per
model instance and examination.

```
cd ~/MCC26/MCC-drivers
./run_oar.sh 'oracle/AirplaneLD-*.out'                      # the warmup shape
TIMEOUT=900 WALLTIME=0:32:0 ./run_oar.sh 'oracle/*-RF.out'   # a real campaign
```

`TIMEOUT` is the per test budget handed to `run_test.pl`, `WALLTIME` the OAR
limit and must exceed it comfortably, `CORES` the cores per job (4, the MCC
count), `HOSTS` an OAR host pattern such as `tall%`, `BK_TOOL` the tool folder.

Two details the script exists for. It submits from a directory named after the
examination suffix, because OAR writes `OAR.<jobid>.stdout` into the directory
the job was submitted from: the results sort themselves into `RC/`, `RF/`,
`LTLF/` and so on. And it sets `RUNATEST=./runatest_cluster.sh`, because a
campaign that spans several examinations has many jobs on the same model
instance running at once, and the default `runatest.sh` would unpack all of
them into one shared `INPUTS/<model>/` that ITS-Tools then writes its unfolded
and reduced nets into. `runatest_cluster.sh` gives each job a `$$` suffixed
copy and removes it at the end.

Submit from the foreground and let it finish; a few hundred `oarsub` calls take
a while, and firing them off in parallel has brought the head down before.

Jobs killed by walltime leave no verdict; resubmit only those by diffing the
oracle list against the tests that reported, as `rerun_oar.sh` does in the
2024 tree:

```
ls oracle/*RF.out | cut -d '/' -f 2 | cut -d '.' -f 1 | sort > oras
grep "Running test" *out | cut -d ':' -f 3 | cut -d '.' -f 2 | sed 's/\s//g' | sort > runs
comm -3 oras runs
```

Watch the queue with `oarstat -u ythierry`; `~/git/ITS-Tools-MCC/monika.sh`
does the same by counting the user's name on the `http://cluster/monika` page.

## Reading the results

The `analysis/` scripts glob `*out`, which matches `OAR.<jobid>.stdout`, so
they are run from inside a per-examination directory and print one CSV to
stdout.

| script | one row per | what it extracts |
| --- | --- | --- |
| `logs2csv.pl` | log | tests started/failed/finished, suite duration, ITS time and memory, verdict counts per technique |
| `itslog2csv.pl` | log | about fifty columns: colored and unfolded net sizes, reduced sizes, counter-example length, invariant counts, traps, **Random/Best/Prob/Parikh walk answers**, every reduction rule, solver attribution |
| `logsITS2csv.pl` | formula | verdict, wall time, ITS time |
| `logs2matrix.pl` | run, on stdin | a `0`/`1`/`.` bitmap of verdicts indexed by formula number, to diff two runs |
| `graphdata.pl`, `graph.sh` | -- | log-log scatter of one column, method against method, through gnuplot |

`itslog2csv.pl` is the one that isolates the PetriSpot contribution: its walk
columns count the properties answered by a walk rather than by SMT or decision
diagrams. It matches ITS-Tools log lines by regex, so the patterns must be
checked against the current output before trusting a campaign.

Our own collector is `Petri/test/mcc/mcclogs2csv.py`: it takes the
per-examination directories and writes one row per log and one row per formula,
carrying the oracle comparison, the failure signature and the PetriSpot walker
census. `Petri/test/mcc/csv/<date>/README.md` reports what a campaign showed.

## Known harness artefact: StateSpace

ITS-Tools reports three of the four `STATE_SPACE` values (`STATES`,
`MAX_TOKEN_IN_PLACE`, `MAX_TOKEN_PER_MARKING`) and not `TRANSITIONS`, the edge
count. Every value it does report matches the oracle, but `run_test.pl` fails
the suite on `3 / 4` results. The exemption for exactly this case is present in
`run_test.pl`, commented out. So an `SS` failure is expected and is not a
regression; check the reported values instead of the exit code.
