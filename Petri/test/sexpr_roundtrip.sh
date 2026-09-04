#!/usr/bin/env bash
# sexpr_roundtrip.sh: check the s-expression property syntax against the MCC
# parser. For each model folder given (model.pnml plus Reachability*.xml):
# MCC XML -> sexpr (names) -> parse -> sexpr must be identical, likewise with
# indices, and the infix dump of the two parses must agree.
#
# Usage: bash Petri/test/sexpr_roundtrip.sh [petri64] bench/models/<model>...
set -u
PETRI=${PETRI:-build/petri64}
if [ -f "${1:-}" ] && [ -x "${1:-}" ] ; then PETRI=$1; shift; fi
LOGS=Petri/test/logs; mkdir -p "$LOGS"
strip() { grep -v '^\[\|^Total runtime\|^;' ; }
fail=0
for dir in "$@"; do
  m=$(basename "$dir")
  for ex in ReachabilityCardinality ReachabilityFireability UpperBounds; do
    xml="$dir/$ex.xml"; [ -f "$xml" ] || continue
    a="$LOGS/$m-$ex.a.sexpr"; b="$LOGS/$m-$ex.b.sexpr"
    ai="$LOGS/$m-$ex.ai.sexpr"; bi="$LOGS/$m-$ex.bi.sexpr"
    "$PETRI" -i "$dir/model.pnml" --props="$xml" --printProps=sexpr | strip > "$a"
    "$PETRI" -i "$dir/model.pnml" --props="$a" --printProps=sexpr | strip > "$b"
    "$PETRI" -i "$dir/model.pnml" --props="$xml" --printProps=sexpr-index | strip > "$ai"
    "$PETRI" -i "$dir/model.pnml" --props="$ai" --printProps=sexpr-index | strip > "$bi"
    "$PETRI" -i "$dir/model.pnml" --props="$xml" --printProps | strip > "$LOGS/$m-$ex.xml.infix"
    "$PETRI" -i "$dir/model.pnml" --props="$a" --printProps | strip > "$LOGS/$m-$ex.sexpr.infix"
    n=$(grep -c '^(' "$a")
    if cmp -s "$a" "$b" && cmp -s "$ai" "$bi" && cmp -s "$LOGS/$m-$ex.xml.infix" "$LOGS/$m-$ex.sexpr.infix" && [ "$n" -gt 0 ]; then
      echo "OK   $m $ex ($n properties)"
    else
      echo "FAIL $m $ex"; fail=1
    fi
  done
done
exit $fail
