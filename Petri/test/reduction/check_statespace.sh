#!/usr/bin/env bash
# check_statespace.sh: one model through `petri64 reduce --goal STATESPACE`, its PNET counted by
# hsc-pn --states, the four values against the deployed StateSpace oracle, 15 s in all.
#
#   bash Petri/test/reduction/check_statespace.sh MODEL [petri64 reduce options...]
#
# MODEL is a folder holding model.pnml (bench/models/X or an INPUTS folder). Prints the
# reduction and counting-record lines, then one line per value: agree, WRONG, or missing
# (a value the consumer left unanswered, as it must when the record does not vouch for it).
# Exit 1 on a WRONG value. Log kept in Petri/test/logs/ss_<model>.log, the export beside it.
set -u
MODEL=$1; shift
ORACLE=${ORACLE:-/data/ythierry/MCC26deploy/MCC-drivers/oracle}
BIN=${BIN:-build/petri64}
HSC=${HSC:-$HOME/git/libHSC/build/tools/hsc-pn}
name=$(basename "$MODEL")
log=Petri/test/logs/ss_${name}.log; out=Petri/test/logs/ss_${name}
rm -rf "$out"
timeout 15s $BIN reduce -i "$MODEL/model.pnml" --goal STATESPACE --output "$out" "$@" > "$log" 2>&1
echo "reduce exit $? ($log)"
grep "^Reduction\|^Counting" "$log"
timeout 15s $HSC --net "$out/model.pnet" --states --totalTime 12 >> "$log" 2>&1
echo "hsc-pn exit $?"
bad=0
for v in STATES TRANSITIONS MAX_TOKEN_IN_PLACE MAX_TOKEN_PER_MARKING; do
  ours=$(awk -v v="$v" '$1=="STATE_SPACE" && $2==v {print $3}' "$log")
  expected=$(awk -v v="$v" '$1=="STATE_SPACE" && $2==v {print $3}' "$ORACLE/$name-SS.out")
  if [ -z "$ours" ]; then s=missing; elif [ "$ours" = "$expected" ]; then s=agree; else s=WRONG; bad=1; fi
  printf '%-8s %-22s %s (oracle %s)\n' "$s" "$v" "${ours:--}" "${expected:-?}"
done
exit $bad
