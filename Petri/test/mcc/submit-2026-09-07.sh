#! /bin/bash
# Campaign of 2026-09-07: the CI product of ITS-Tools bae71d95 (Effort contract, spotutil,
# Spot 2.16) with the petri64 of PetriSpot 3f29014 (step cap of a focused round).
# One run_oar.sh after another, never in parallel; the result folders of the previous
# campaign were archived to /data/ythierry/MCC26archive/2026-09-06c and removed first.
# Runs on the cluster head from ~/MCC26/MCC-drivers, detached:
#   setsid nohup ./submit-2026-09-07.sh > submit-2026-09-07.log 2>&1 &
cd ~/MCC26/MCC-drivers
EXAMS="RD QLA LTLC LTLF"
for ex in $EXAMS ; do
	if [ -d $ex ] ; then mv $ex $ex-before-2026-09-07 ; fi
done
for ex in $EXAMS ; do
	date
	TIMEOUT=1800 WALLTIME=0:45:00 HOSTS=tall% CORES=4 RUNATEST=./runatest_cluster.sh ./run_oar.sh "oracle/*-$ex.out"
done
date
echo SUBMISSION DONE
