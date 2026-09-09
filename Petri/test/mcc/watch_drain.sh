#!/usr/bin/env bash
# Watch our jobs on the cluster until none is queued or running: one
# cluster_status.sh every INTERVAL seconds (600), each line appended to LOG,
# exit 0 when the queue has drained, exit 1 after MAXLOOKS looks (600).
#
#   watch_drain.sh LOG [INTERVAL] [MAXLOOKS] [EXAM...]
set -u
LOG=${1:?log file}; INTERVAL=${2:-600}; MAXLOOKS=${3:-600}; shift; shift 2>/dev/null; shift 2>/dev/null
DIR=$(cd "$(dirname "$0")" && pwd)
for i in $(seq 1 "$MAXLOOKS"); do
  out=$(bash "$DIR/cluster_status.sh" "$@" 2>&1)
  echo "$out" >> "$LOG"
  jobs=$(echo "$out" | grep '^jobs:' | sed 's/jobs: *//')
  if [ -z "$(echo "$jobs" | tr -d ' ')" ]; then
    echo "drained after $i looks ($(date))" | tee -a "$LOG"
    exit 0
  fi
  sleep "$INTERVAL"
done
echo "still running after $MAXLOOKS looks" | tee -a "$LOG"
exit 1
