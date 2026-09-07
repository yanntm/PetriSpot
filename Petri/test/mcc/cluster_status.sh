#!/usr/bin/env bash
# cluster_status.sh: how far a campaign is on the cluster, in one line per examination.
#
#   bash Petri/test/mcc/cluster_status.sh [EXAM...]     # default: every result folder present
#
# Prints, per examination folder of ~/MCC26/MCC-drivers on the cluster, the logs written
# (one OAR.<id>.stdout per finished job) and how many of them carry the harness trailer
# (a job killed at the wall leaves a log without it); then the user's OAR jobs by state
# (R running, W waiting) and the cluster's time. Read-only.
set -u
HOST=${CLUSTER_HOST:-cluster.lip6.fr}
TREE=${CLUSTER_TREE:-MCC26/MCC-drivers}
EXAMS="$*"
ssh -o BatchMode=yes -o ConnectTimeout=20 "$HOST" "cd $TREE || exit 1
exams=\"$EXAMS\"
if [ -z \"\$exams\" ]; then exams=\$(ls -d RC RF RD UB L QL SM OS SS LTLC LTLF CTLC CTLF QLA SMA UBA 2>/dev/null | tr '\n' ' '); fi
for d in \$exams; do
  [ -d \$d ] || { echo \"\$d: no folder\"; continue; }
  logs=\$(ls \$d | grep -c 'stdout\$')
  done=\$(grep -l 'Test :' \$d/*.stdout 2>/dev/null | wc -l)
  echo \"\$d: \$logs logs, \$done with a trailer\"
done
echo \"jobs: \$(oarstat -u 2>/dev/null | tail -n +3 | awk '{print \$2}' | sort | uniq -c | awk '{printf \"%s %s  \", \$1, \$2}')\"
date"
