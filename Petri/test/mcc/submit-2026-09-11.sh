#! /bin/bash
# Campaign of 2026-09-11: a fresher ITS-Tools baseline on the two reachability
# examinations. The CI product 202609101624 (its bundled engines: the hsc-pn of
# the libHSC CI, the petri64 of the PetriSpot CI), launched through a GraalVM
# native image built locally for x86-64-v2 (CLUSTER.md §1, "Which launcher a
# campaign gets"): the CI's own image is x86-64-v3 and answers nothing on the
# Westmere nodes, and this campaign runs on small% because tall% is saturated.
# 1800 s a run, 6 cores: OAR caps a job's memory at the node's RAM over its
# declared CPUs times the cores asked, and small% is 64 GB over 24, so 6 cores
# is the 16 GB of the MCC rule.
# TAG=itstools writes RC.itstools / RF.itstools, beside the older tag-less
# folders of the previous campaigns; RUNATEST=./runatest_cluster.sh gives each
# job its own copy of the model, as the structural reductions write there.
# One run_oar.sh after another, never in parallel.
# Runs on the cluster head from ~/MCC26/MCC-drivers, detached:
#   setsid nohup ./submit-2026-09-11.sh > submit-2026-09-11.log 2>&1 &
cd ~/MCC26/MCC-drivers
EXAMS="RC RF"
for ex in $EXAMS ; do
	date
	TIMEOUT=1800 WALLTIME=0:45:0 HOSTS=small% CORES=6 TAG=itstools \
	  RUNATEST=./runatest_cluster.sh BK_TOOL=itstools ./run_oar.sh "oracle/*-$ex.out"
done
date
echo SUBMISSION DONE
