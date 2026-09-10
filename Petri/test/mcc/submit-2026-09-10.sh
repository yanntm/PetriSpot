#! /bin/bash
# Campaign of 2026-09-10: libHSC alone (the hsc tool folder, hsc-pn from libHSC
# 142fd4b: a stopped closure decides nothing, a partial reachable set answers only
# what stands, a deadline is never an error) on the CTL examinations of the PT
# nets, the only ones the tool declares -- the first CTL run of hsc-pn since the
# order sweep found its wrong verdicts. On small% at 6 cores: the nodes are
# x86-64-v2 (no AVX, no native image needed here) and OAR caps a job's memory at
# the node's RAM per core times the cores asked, 6 cores being 16 GB on small.
# 600 s like the hsc600 campaign, the tag hsc as it writes CTLC.hsc / CTLF.hsc.
# One run_oar.sh after another, never in parallel; the warmup folders of the
# same day are moved aside first.
# Runs on the cluster head from ~/MCC26/MCC-drivers, detached:
#   setsid nohup ./submit-2026-09-10.sh > submit-2026-09-10.log 2>&1 &
cd ~/MCC26/MCC-drivers
EXAMS="CTLC CTLF"
for ex in $EXAMS ; do
	if [ -d $ex.hsc ] ; then mv $ex.hsc $ex.hsc-warmup-2026-09-10 ; fi
done
for ex in $EXAMS ; do
	date
	TIMEOUT=600 WALLTIME=0:15:0 HOSTS=small% CORES=6 TAG=hsc BK_TOOL=hsc ./run_oar.sh "oracle/*-PT-*-$ex.out"
done
date
echo SUBMISSION DONE
