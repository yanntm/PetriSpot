#!/usr/bin/env bash
# ctl_oracle.sh: the CTL checker against the contest oracles, several seeds.
# For each model folder given, CTLCardinality and CTLFireability are run with
# each seed under --totalTime and every FORMULA line is compared with
# ~/git/MCC-drivers/oracle/<model>-CTLC.out / -CTLF.out. Prints one line per
# run (ok / unknown / WRONG counts) and every WRONG verdict; logs in
# Petri/test/logs/<model>-<ex>-seed<k>.log.
#
# Usage: bash Petri/test/ctl_oracle.sh [-t seconds] [-s "1 2 3"] bench/models/<model>...
set -u
PETRI=${PETRI:-build/petri64}
ORACLE=${ORACLE:-$HOME/git/MCC-drivers/oracle}
TOTAL=10; SEEDS="1 2 3"
while getopts "t:s:" opt; do case $opt in t) TOTAL=$OPTARG;; s) SEEDS=$OPTARG;; esac; done
shift $((OPTIND-1))
LOGS=Petri/test/logs; mkdir -p "$LOGS"
wrong=0
for dir in "$@"; do
  m=$(basename "$dir")
  for ex in CTLCardinality:CTLC CTLFireability:CTLF; do
    x=${ex%%:*}; short=${ex##*:}
    xml="$dir/$x.xml"; [ -f "$xml" ] || continue
    orc="$ORACLE/$m-$short.out"; [ -f "$orc" ] || { echo "no oracle for $m $short"; continue; }
    for seed in $SEEDS; do
      log="$LOGS/$m-$x-seed$seed.log"
      "$PETRI" -i "$dir/model.pnml" --props="$xml" --totalTime=$TOTAL --seed=$seed --printUnknown -q --trace > "$log" 2>&1
      ok=0; unk=0; bad=0
      while read -r name verdict; do
        ours=$(grep -E "^(FORMULA|UNKNOWN) $name( |$)" "$log" | head -1 | awk '{print $1=="UNKNOWN"?"?":$3}')
        if [ "$ours" = "$verdict" ]; then ok=$((ok+1));
        elif [ -z "$ours" ] || [ "$ours" = "?" ]; then unk=$((unk+1));
        else bad=$((bad+1)); echo "WRONG $m $x seed $seed: $name oracle $verdict ours $ours"; fi
      done < <(grep "^FORMULA" "$orc" | awk '{print $2, $3}')
      echo "$m $x seed $seed: $ok ok, $unk unknown, $bad wrong ($(grep '^Total' "$log"))"
      [ $bad -gt 0 ] && wrong=1
    done
  done
done
exit $wrong
