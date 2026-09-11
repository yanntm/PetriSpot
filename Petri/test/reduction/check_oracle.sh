#!/usr/bin/env bash
# check_oracle.sh: one model, one examination, petri64 against the deployed oracle, 15 s.
#
#   bash Petri/test/reduction/check_oracle.sh MODEL EXAM [petri64 options...]
#
# MODEL is a folder holding model.pnml and the MCC formula files (bench/models/X or an
# INPUTS folder), EXAM one of RC RF CTLC CTLF UB RD. Runs build/petri64 on the formula
# file with the options given (--reduce, -q, budgets...) under a 15 s cap, then prints the
# reduction and initial-state lines, and per FORMULA line whether the deployed oracle
# agrees. Exit 1 on any disagreement. Log kept in Petri/test/logs/check_<model>_<exam>.log.
set -u
MODEL=$1; EXAM=$2; shift 2
ORACLE=${ORACLE:-/data/ythierry/MCC26deploy/MCC-drivers/oracle}
BIN=${BIN:-build/petri64}
name=$(basename "$MODEL")
case $EXAM in
  RC) file=ReachabilityCardinality.xml;; RF) file=ReachabilityFireability.xml;;
  CTLC) file=CTLCardinality.xml;; CTLF) file=CTLFireability.xml;; UB) file=UpperBounds.xml;;
  RD) file=;; *) echo "unknown examination $EXAM" >&2; exit 2;;
esac
log=Petri/test/logs/check_${name}_${EXAM}.log
if [ -z "$file" ]; then timeout 15s $BIN -i "$MODEL/model.pnml" --findDeadlock "$@" > "$log" 2>&1
else timeout 15s $BIN -i "$MODEL/model.pnml" --props="$MODEL/$file" "$@" > "$log" 2>&1; fi
echo "exit $? ($log)"
grep -i "^Reduction\|Initial state" "$log"
grep "^FORMULA" "$log" | awk '{print $2, $3, $6, $7}' | sort > "$log.ours"
awk '/^FORMULA/{print $2, $3}' "$ORACLE/$name-$EXAM.out" | sort > "$log.oracle"
bad=0
while read -r id verdict t1 t2; do
  expected=$(awk -v id="$id" '$1==id{print $2}' "$log.oracle")
  if [ "$expected" = "$verdict" ]; then s=agree; elif [ "$expected" = "?" ] || [ -z "$expected" ]; then s=unverified; else s=WRONG; bad=1; fi
  printf '%-8s %-58s %-6s %s %s\n' "$s" "$id" "$verdict" "$t1" "$t2"
done < "$log.ours"
echo "answered $(wc -l < "$log.ours") / $(wc -l < "$log.oracle")"
rm -f "$log.ours" "$log.oracle"
exit $bad
