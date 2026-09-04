#!/usr/bin/env bash
# pnet_roundtrip.sh: the PNET path must behave as the PNML path. For each
# model folder given: export model.pnml as PNET, convert the MCC properties to
# index-based s-expressions, then compare (fixed seed, random strategy) the
# FORMULA lines and step counts of PNML+XML against PNET+sexpr, and the
# invariant counts of both inputs.
#
# Walks are step-bound (--walkSteps) so that they are reproducible; keep to
# small models in the development loop.
#
# Usage: bash Petri/test/pnet_roundtrip.sh [petri64] bench/models/<model>...
set -u
PETRI=${PETRI:-build/petri64}
if [ -f "${1:-}" ] && [ -x "${1:-}" ] ; then PETRI=$1; shift; fi
LOGS=Petri/test/logs; mkdir -p "$LOGS"
walk() { grep '^FORMULA\|^Thread' | sed 's/, [0-9]* ms (.*//' ; }
fail=0
for dir in "$@"; do
  m=$(basename "$dir"); pnet="$LOGS/$m.pnet"
  "$PETRI" -i "$dir/model.pnml" --exportNet="$pnet" -q > /dev/null || { echo "FAIL $m export"; fail=1; continue; }
  a=$("$PETRI" -i "$dir/model.pnml" --Pflows --Tflows -q | grep '^Computed' | sed 's/ in .*//')
  b=$("$PETRI" --net="$pnet" --Pflows --Tflows -q | grep '^Computed' | sed 's/ in .*//')
  if [ "$a" == "$b" ] && [ -n "$a" ]; then echo "OK   $m invariants"; else echo "FAIL $m invariants"; fail=1; fi
  for ex in ReachabilityCardinality ReachabilityFireability; do
    xml="$dir/$ex.xml"; [ -f "$xml" ] || continue
    sx="$LOGS/$m-$ex.pnet.sexpr"
    "$PETRI" -i "$dir/model.pnml" --props="$xml" --printProps=sexpr-index | grep '^(' > "$sx"
    "$PETRI" -i "$dir/model.pnml" --props="$xml" --seed=7 --strategy=random -t 20 --walkSteps=20000 -q | walk > "$LOGS/$m-$ex.pnml.walk"
    "$PETRI" --net="$pnet" --props="$sx" --seed=7 --strategy=random -t 20 --walkSteps=20000 -q | walk > "$LOGS/$m-$ex.pnet.walk"
    n=$(grep -c '^FORMULA' "$LOGS/$m-$ex.pnml.walk")
    if cmp -s "$LOGS/$m-$ex.pnml.walk" "$LOGS/$m-$ex.pnet.walk"; then echo "OK   $m $ex ($n verdicts)"; else echo "FAIL $m $ex"; fail=1; fi
  done
done
exit $fail
