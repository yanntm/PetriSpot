#! /bin/bash
# Rate of the wrong verdict on StigmergyCommit-PT-02b-LTLCardinality-03 (oracle FALSE).
# The property is closed by whichever engine gets there first, so one run says nothing.
#   bash stigmergy-repeat.sh <n> <product folder> [more product folders...]
set -u
N=${1:?number of runs}; shift
M=/home/ythierry/git/PetriSpot/bench/models/Stigmergy-f03
G=/data/ythierry/MCC26deploy/MCC-drivers/itstools/greatspn/
L=/home/ythierry/git/PetriSpot/Petri/test/logs
for P in "$@" ; do
	echo "===== $(ls -d $P/plugins/fr.lip6.move.gal.application.pnmcc_*.jar | sed 's/.*pnmcc_//;s/\.jar//')  ($P)"
	for i in $(seq 1 $N) ; do
		o=$L/stig-$(basename $P)-$i.log
		( cd "$P" && ITSTOOLS="$P" ./its-tools-flat.sh -pnfolder $M -examination LTLCardinality \
			-its -ltsmin -greatspnpath $G -order META -manyOrder -smt -timeout 1800 ) > $o 2>&1
		printf "  %2d  %s\n" $i "$(grep -m1 '^FORMULA' $o | cut -d' ' -f3-)"
	done
done
