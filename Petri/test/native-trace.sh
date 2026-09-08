#! /bin/bash
# Record GraalVM reachability metadata for the examinations the ITS-Tools native image
# config never covered. The agent merges into ~/git/ITStools/ITS-commandline/native/config/;
# the image is then rebuilt from it by the ITS-Tools CI.
#
# The recorded config came from AirplaneLD PT-0010 (OneSafe, deadlock, LTLC, UB) and
# COL-0010 (LTLF, CTLF, RC). StateSpace was never traced and its decomposition path
# (MultiOrderRunner.runMultiITS -> CompositeBuilder) misses fr.lip6.move.gal.InstanceDecl[].
#
# One run per pair, each allowed to fail: trace-agent.sh has set -e and would stop the
# whole sweep on the first non-zero exit.
set -u
export GRAALVM=${GRAALVM:-/data/ythierry/graal/graalvm-jdk-25.0.4+7.1}
NATIVE=${NATIVE:-$HOME/git/ITStools/ITS-commandline/native}
PROD=${PROD:-/data/ythierry/MCC26deploy/MCC-drivers/itstools/itstools}
IN=${IN:-/data/ythierry/MCC26deploy/MCC-drivers/INPUTS}
PT="$IN/AirplaneLD-PT-0010"
COL="$IN/AirplaneLD-COL-0010"
# AirplaneLD-COL declares no product sort, so the symmetric net terms of a
# coloured model with one were never recorded: BART-COL died on
# fr.lip6.move.pnml.symmetricnet.terms.Sort[] in the campaign of 2026-09-08.
BART="$IN/BART-COL-002"
# Neither Airplane nor BART carries an XML comment, so the Axiom node the PNML
# framework builds for one was never instantiated: the COL models the GreatSPN
# Editor exported (Sudoku, LastZero, FileSystem, UtilityControlRoom, VehicularWifi)
# all died at 11 ms on org.apache.axiom.om.impl.llom.OMCommentImpl in the campaign
# of 2026-09-08, 58 instances answering nothing.
SUDOKU="$IN/Sudoku-COL-AN01"
run() {
	echo "===== trace $2 on $(basename $1)"
	"$NATIVE/trace-agent.sh" "$PROD" "$1" "$2" > /dev/null 2>&1 || echo "  (exit $?, traced anyway)"
}
for e in StateSpace Liveness CTLCardinality QuasiLiveness StableMarking ReachabilityFireability ; do
	run "$PT" "$e"
done
for e in StateSpace Liveness CTLCardinality ; do
	run "$COL" "$e"
done
for e in CTLCardinality CTLFireability ReachabilityCardinality LTLFireability ; do
	run "$BART" "$e"
done
for e in CTLCardinality CTLFireability ; do
	run "$SUDOKU" "$e"
done
echo "===== done; config entries now:"
python3 -c "import json;d=json.load(open('$NATIVE/config/reachability-metadata.json'));print({k:len(v) for k,v in d.items()})"
