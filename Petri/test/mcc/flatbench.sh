#! /bin/bash
# Startup comparison on one node: the Eclipse launcher, the flat classpath and the native image, OneSafe on a small net.
cd ~/MCC26/flat-test
I=~/MCC26/MCC-drivers/INPUTS/AirplaneLD-PT-0010
G=~/MCC26/MCC-drivers/itstools/greatspn/
uname -n
for L in eclipse flat native eclipse flat native eclipse flat native ; do
	if [ $L = flat ] ; then CMD=./its-tools-flat.sh ; elif [ $L = native ] ; then CMD="./its-tools-native -Dfr.lip6.binaries.root=$PWD/plugins" ; else CMD=./its-tools ; fi
	S=$(date +%s%N)
	$CMD -pnfolder $I -examination OneSafe -its -ltsmin -greatspnpath $G -order META -manyOrder -smt -timeout 300 > run-$L.log 2>&1
	E=$(date +%s%N)
	echo "$L: wall $(( (E-S)/1000000 )) ms, internal $(grep -o 'Total runtime [0-9]* ms' run-$L.log | tail -1), $(grep -c FORMULA run-$L.log) FORMULA, $(grep -ci 'exception\|error has occurred' run-$L.log) errors"
done
S=$(date +%s%N); java -Xss128m -Xmx16384m -version > /dev/null 2>&1; E=$(date +%s%N); echo "bare JVM: $(( (E-S)/1000000 )) ms"
