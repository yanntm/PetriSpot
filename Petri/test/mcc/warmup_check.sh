#!/usr/bin/env bash
# Read a warmup: rsync the named result folders of the cluster tree down (never
# --delete) into $DEST and print, per folder, the logs, the trailers, the runs
# with a regression line (a wrong or missing value against the oracle), the
# runs with a failure signature, and the FORMULA / STATE_SPACE lines counted.
#
#   warmup_check.sh DEST FOLDER...      e.g. warmup_check.sh /data/ythierry/MCC26run/warmup-2026-09-10 RC RF CTLC.hsc RC.hscapprox
set -u
DEST=${1:?destination}; shift
HOST=${CLUSTER_HOST:-cluster.lip6.fr}
TREE=${CLUSTER_TREE:-MCC26/MCC-drivers}
mkdir -p "$DEST"
for d in "$@"; do
  rsync -rlptD --no-g -q "$HOST:$TREE/$d/" "$DEST/$d/" 2>/dev/null || { echo "$d: no folder on the cluster"; continue; }
  logs=$(ls "$DEST/$d" 2>/dev/null | grep -c 'stdout$')
  trailers=$(grep -l 'Test :' "$DEST/$d"/*.stdout 2>/dev/null | wc -l)
  regress=$(grep -l 'regression detected\|testFailed' "$DEST/$d"/*.stdout 2>/dev/null | wc -l)
  fails=$(grep -l 'Segmentation\|bad_alloc\|Exception\|CANNOT_COMPUTE\|No such file' "$DEST/$d"/*.stdout 2>/dev/null | wc -l)
  answers=$(cat "$DEST/$d"/*.stdout 2>/dev/null | grep -c '^FORMULA \|^STATE_SPACE ')
  echo "$d: $logs logs, $trailers trailers, $regress with a regression line, $fails with a failure signature, $answers answers"
  grep -h 'testFailed' "$DEST/$d"/*.stdout 2>/dev/null | cut -c1-200 | head -5
done
