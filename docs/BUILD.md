# Guide: building and running locally

## PetriSpot

```
cmake -S Petri -B build && cmake --build build --target petri64     # the development loop, about a minute after a header change
./buildPetriSpot.sh                                                  # what the CI runs: static libexpat, stripped petri32/64/128 in website/
```

Binaries: `build/petri32|64|128` (integer width) and `build/kersconv`.
Other trees in use: `build-prof/` (profiling flags), `build-tsan/` (thread
sanitizer). A header-only project: touching `Petri/src/**` rebuilds the one
translation unit `Petri.cpp`.

Quick checks after a change:

```
./build/petri64 -i Petri/examples/Airplane.pnml --props Petri/test/props/Airplane.sexpr -q -t 3        # reachability forms
./build/petri64 -i Petri/examples/Airplane.pnml --props Petri/test/props/Airplane-ctl.sexpr --trace -t 3   # CTL forms, with evidence
bash Petri/test/sexpr_roundtrip.sh bench/models/AirplaneLD-PT-0010                                    # MCC XML <-> s-expressions
bash Petri/test/ctl_oracle.sh -t 8 -s "1 2" bench/models/AirplaneLD-PT-0010                            # CTL verdicts against the oracle
```

Development models live in `bench/models/<model>/` (git-ignored): extract
`~/git/pnmcc-models-2026/website/INPUTS/<model>.tgz` there. Keep every run
under about 15 s and redirect long output to `Petri/test/logs/`.

## ITS-Tools

```
cd ~/git/ITStools/fr.lip6.move.gal.parent && mvn -o install -DskipTests > ~/git/PetriSpot/Petri/test/logs/itstools-mvn-<tag>.log 2>&1   # read, then rm
```

About 1:30 min offline once the Tycho cache is warm (`-o`); the product
tarball lands in
`ITS-commandline/fr.lip6.move.gal.itscl.product/target/products/fr.lip6.move.gal.itscl.product-linux.gtk.x86_64.tar.gz`.
The build fetches `petri64` from `Inv-Linux`, so a local product carries the
*published* PetriSpot, not the working tree. To run a local `petri64` inside
a local product:

```
P=/data/ythierry/MCC26deploy/products/itstools-local && rm -rf $P && mkdir -p $P   # a local product lives under MCC26deploy/products/, never at the top of /data/ythierry
tar -xzf ~/git/ITStools/ITS-commandline/fr.lip6.move.gal.itscl.product/target/products/fr.lip6.move.gal.itscl.product-linux.gtk.x86_64.tar.gz -C $P
cp ~/git/PetriSpot/build/petri64 $P/plugins/fr.lip6.petrispot.binaries_*/bin/petri64
```

(or the JVM property `-Dpetrispot.bin=<path>`, read by `PetriSpotWalker`).
The tree `~/git/ITStools` carries the user's own uncommitted edits: add only
the files you changed, never `git add -A`.

## Running ITS-Tools on a model folder

```
cd bench/models/AirplaneLD-PT-0010
$P/its-tools -pnfolder . -examination CTLFireability -its -smt -timeout 120 > log 2>&1
$P/its-tools-flat.sh -pnfolder . -examination LTLCardinality -its -timeout 120     # no Equinox, faster start
```

The examination names are the MCC ones (`ReachabilityCardinality`,
`Liveness`, `CTLFireability`...). `-its` enables the decision diagram
engines, `-smt` the SMT presolving, `-ltsmin -greatspnpath <dir>` the other
engines the contest driver adds (`ITS-Tools-MCC/runeclipse.sh` has the full
line: `-its -ltsmin -greatspnpath $BINDIR/greatspn/ -order META -manyOrder
-smt -timeout $BK_TIME_CONFINEMENT`). What the run printed about PetriSpot:

```
grep -E "^FORMULA|PetriSpot walker|beside the decision diagrams|Total runtime" log
```

`PetriSpotWalker.DEBUG` (1 keeps the exchanged files, 2 echoes the binary's
output) is the switch when a call has to be seen.

## Through the harness

`docs/CLUSTER.md` section 6: `run_test.pl` against an oracle file, the same
log format as the cluster.
