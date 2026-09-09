#! /bin/bash
# Campaign of 2026-09-10, three legs, one examination per batch, the queue
# drained between batches (never more than one examination's jobs queued):
#   1. itstools  RC, RF        1800 s  the unvalidated ITS-Tools / PetriSpot changes,
#                                      on the deployed product (no libHSC in that chain)
#   2. hsc       CTLC, CTLF    1800 s  libHSC's CTL checker alone, four shapes, P/T only
#   3. hscapprox RC, RF         300 s  libHSC with the over-approximation first, P/T only
# libHSC binaries: the HSC-Linux build of 24c28ad (hsc/bin, shared by hscapprox).
# Runs on the cluster head from ~/MCC26/MCC-drivers, detached:
#   setsid nohup ./submit-2026-09-10.sh > submit-2026-09-10.log 2>&1 &
# Warmups (AirplaneLD-PT-0010 and a handful of instances per leg) were run
# and read before this script was started.
cd ~/MCC26/MCC-drivers || exit 1

drained() {  # wait until none of our jobs is queued or running; 10 minutes between looks
	while true ; do
		n=$(oarstat -u 2>/dev/null | tail -n +3 | grep -c .)
		if [ "$n" -eq 0 ] ; then return 0 ; fi
		echo "$(date) $n jobs still in the queue"
		sleep 600
	done
}

batch() {  # tool exam pattern timeout walltime [tag]
	local tool=$1 ex=$2 pat=$3 to=$4 wt=$5 tag=${6:-}
	date
	echo "== $tool $ex ($pat) timeout $to tag '$tag'"
	TAG=$tag TIMEOUT=$to WALLTIME=$wt HOSTS=tall% CORES=4 RUNATEST=./runatest_cluster.sh BK_TOOL=$tool ./run_oar.sh "$pat"
	sleep 120
	drained
}

for ex in RC RF ; do
	if [ -d $ex ] ; then mv $ex $ex-before-2026-09-10 ; fi
	batch itstools $ex "oracle/*-$ex.out" 1800 0:45:00
done
for ex in CTLC CTLF ; do
	if [ -d $ex.hsc ] ; then mv $ex.hsc $ex.hsc-before-2026-09-10 ; fi
	batch hsc $ex "oracle/*-PT-*-$ex.out" 1800 0:45:00 hsc
done
for ex in RC RF ; do
	if [ -d $ex.hscapprox ] ; then mv $ex.hscapprox $ex.hscapprox-before-2026-09-10 ; fi
	batch hscapprox $ex "oracle/*-PT-*-$ex.out" 300 0:10:00 hscapprox
done
date
echo SUBMISSION DONE
