#!/usr/bin/env bash
# PT-only warmups and serial bundled campaigns. Deploy unchanged before use.
# warmup: 3 StateSpace portfolio jobs, then 2 models x SS/CTLC/CTLF sweep jobs.
# full: requires the locally reviewed warmup marker; never run concurrently.
set -euo pipefail
MODE=${1:?warmup or full}
TREE=$HOME/MCC26/MCC-drivers
SWEEP=$HOME/MCC26/hsc-sweep-reduce-20260912
TAG=reduce20260912
cd "$TREE"
exec 9>"submit-$TAG.lock"
flock -n 9 || { echo 'Another submission is active'; exit 1; }
queue_empty() { [ -z "$(oarstat -u)" ]; }
wait_capacity() {
  while [ "$(oarstat -u | awk '$1 ~ /^[0-9]+$/ {n++} END {print n+0}')" -gt 3000 ]; do
    sleep 300
  done
}
validate() { [[ "$1" == *-PT-* && "$1" != *-COL-* && "$1" =~ ^[A-Za-z0-9_-]+$ ]]; }
submit_ss() {
  local model=$1 budget=$2 wall=$3 tag=$4
  validate "$model" || { echo "Rejected non-PT model: $model"; exit 2; }
  test -f "oracle/$model-SS.out"
  test -f "INPUTS/$model.tgz"
  mkdir -p "SS.$tag"
  (
    cd "SS.$tag"
    oarsub -n "SS-$tag" -l "/nodes=1/core=6,walltime=$wall" -p "host like 'small%'" \
      "cd $TREE && HSC_REDUCE=1 BK_TOOL=hsc RUNATEST=./runatest_cluster.sh ./run_test.pl oracle/$model-SS.out -t $budget"
  )
  echo "submitted SS $model"
}
case "$MODE" in
  warmup)
    queue_empty || { echo 'Queue not empty; no warmup submitted'; exit 1; }
    test ! -e "warmup-$TAG.started"
    date > "warmup-$TAG.started"
    for model in AirplaneLD-PT-0010 TokenRing-PT-005 DoubleExponent-PT-003; do
      submit_ss "$model" 60 0:03:00 "warmup-$TAG"
    done
    cd "$SWEEP"
    for ex in SS CTLC CTLF; do
      CORES=6 SWEEP_REDUCE=1 ./sweep_oar.sh -x "$ex" -t 300 -w 1:34:00 -H 'small%' -o "warmup-$TAG" models_warmup.txt
    done
    ;;
  full)
    test -f "warmup-$TAG.reviewed" || { echo 'Warmups have not been reviewed'; exit 1; }
    queue_empty || { echo 'Queue not empty; no campaign submitted'; exit 1; }
    test ! -e "full-$TAG.started"
    date > "full-$TAG.started"
    while read -r model rest; do
      [ -n "$model" ] || continue
      submit_ss "$model" 900 0:17:00 "$TAG"
    done < "$SWEEP/models_pt_all.txt"
    echo 'StateSpace submitted'
    cd "$SWEEP"
    for ex in SS CTLC CTLF; do
      wait_capacity
      CORES=6 SWEEP_REDUCE=1 ./sweep_oar.sh -x "$ex" -t 300 -w 1:34:00 -H 'small%' -o "$TAG" models_pt_all.txt
    done
    ;;
  *) echo "Unknown mode $MODE"; exit 2;;
esac
date
echo 'SUBMISSION DONE'
