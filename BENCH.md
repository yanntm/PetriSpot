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
and `TEDD2026` verdicts), so `install_oracle.sh` is not needed.

`INPUTS/` holds 1951 archives, the whole 2026 set minus `StigmergyCommit-PT-11b`
and `TokenRing-PT-050`, which exceed the GitHub pages file size limit and are
therefore absent from `pnmcc-models-2026` itself.

## Rsync to the cluster

```
rsync -rlptD --no-g --chmod=Dg+s --delete /data/ythierry/MCC26deploy/ cluster.lip6.fr:MCC26/
```

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

One `oarsub` per oracle file, that is per model instance and examination, on
one node and four cores, with a walltime of about twice the `-t` budget. The
submission happens from a per-examination directory because OAR writes
`OAR.<jobid>.stdout` into the directory the job was submitted from, which is
how the results sort themselves:

```
for i in oracle/*-RF.out ; do
	DIR=$(echo $i | perl -pe 's/.*\-(\w+)\.out/\1/g')
	mkdir -p $DIR ; cd $DIR
	oarsub -l "/nodes=1/core=4,walltime=0:32:0" \
	  "cd ~/MCC26/MCC-drivers/ && BK_TOOL=itstools ./run_test.pl $i -t 900 ; exit"
	cd ..
done
```

Jobs killed by walltime leave no verdict; resubmit only those by diffing the
oracle list against the tests that reported, as `rerun_oar.sh` does in the
2024 tree:

```
ls oracle/*RF.out | cut -d '/' -f 2 | cut -d '.' -f 1 | sort > oras
grep "Running test" *out | cut -d ':' -f 3 | cut -d '.' -f 2 | sed 's/\s//g' | sort > runs
comm -3 oras runs
```

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

## Known harness artefact: StateSpace

ITS-Tools reports three of the four `STATE_SPACE` values (`STATES`,
`MAX_TOKEN_IN_PLACE`, `MAX_TOKEN_PER_MARKING`) and not `TRANSITIONS`, the edge
count. Every value it does report matches the oracle, but `run_test.pl` fails
the suite on `3 / 4` results. The exemption for exactly this case is present in
`run_test.pl`, commented out. So an `SS` failure is expected and is not a
regression; check the reported values instead of the exit code.
