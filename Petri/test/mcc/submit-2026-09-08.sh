#! /bin/bash
# Campaign of 2026-09-08: the first CTL examinations at scale. The CI product of
# ITS-Tools 202609071637 (the explicit CTL checker plugged beside its-ctl) with the
# petri64 of PetriSpot 16cb816 (--ctlSteps), launched through the GraalVM native
# image: runeclipse.sh execs its-tools-native when it is present, and the image
# needs AVX2, so tall% only (on small% it refuses and answers nothing).
# Liveness follows the two CTL examinations: it takes the same path, stated as
# AG EF fireable per transition.
# One run_oar.sh after another, never in parallel; the folders of the previous
# campaign were archived to /data/ythierry/MCC26archive/2026-09-07.
# Runs on the cluster head from ~/MCC26/MCC-drivers, detached:
#   setsid nohup ./submit-2026-09-08.sh > submit-2026-09-08.log 2>&1 &
cd ~/MCC26/MCC-drivers
EXAMS="CTLC CTLF L"
for ex in $EXAMS ; do
	if [ -d $ex ] ; then mv $ex $ex-before-2026-09-08 ; fi
done
for ex in $EXAMS ; do
	date
	TIMEOUT=1800 WALLTIME=0:45:00 HOSTS=tall% CORES=4 RUNATEST=./runatest_cluster.sh ./run_oar.sh "oracle/*-$ex.out"
done
date
echo SUBMISSION DONE
