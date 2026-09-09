#!/usr/bin/env bash
# reap_reorg.sh: the one-time move of every collected campaign into the single
# log root, /data/ythierry/MCC26logs/<tool>/<build>/<EXAM>[.tag]/ (CLUSTER.md §4).
# A campaign is labelled by the product build id read from one of its logs
# (the ITS-Tools plugin timestamp, 1.0.0.<id>, or the native image's); a run
# that mixed builds keeps the one grabbed. Two campaigns of one build on the
# same examination: the second folder gets the source name as a tag.
set -eu
R=/data/ythierry/MCC26logs; mkdir -p "$R"
place() {  # src tool build [exam-name-for-loose-folder]
  local src=$1 tool=$2 build=$3 loose=${4:-}; local dst="$R/$tool/$build"; mkdir -p "$dst"
  local tag; tag=$(basename "$src")
  if [ -n "$loose" ]; then
    local d="$dst/$loose"; [ -e "$d" ] && d="$dst/$loose.$tag"; mv "$src" "$d"; echo "  $src -> $d"; return
  fi
  for e in "$src"/* "$src"/.[!.]*; do [ -e "$e" ] || continue; local b; b=$(basename "$e")
    if [ -d "$e" ] && [ -n "$(find "$e" -maxdepth 1 -name 'OAR.*.stdout' -print -quit)" ]; then
      local d="$dst/$b"; [ -e "$d" ] && d="$dst/$b.$tag"; mv "$e" "$d"; echo "  $e -> $d"
    elif [ "$b" = csv ]; then mv "$e" "$dst/csv.$tag"; echo "  $e -> $dst/csv.$tag"
    else mkdir -p "$dst/_$tag"; mv "$e" "$dst/_$tag/"; fi
  done
  rmdir "$src" && echo "  removed empty $src"
}
A=/data/ythierry/MCC26archive; U=/data/ythierry/MCC26run
place $A/2026-09-06-baseline itstools 2026090600
place $U/2026-09-06 itstools 2026090600
place $A/2026-09-06c itstools 2026090616
place $U/2026-09-06c itstools 2026090616
place $U/L itstools 2026090513 L
place $U/QL itstools 2026090513 QL
place $U/UB itstools 2026090513 UB
place $A/2026-09-05 itstools 2026090513
place $U/SM itstools 202609051349 SM
place $A/2026-09-06a itstools 202609051349
place $U/RD itstools 202609052009 RD
place $A/RDvar itstools 202609052009 RD.var
place $A/2026-09-06-precampaign itstools 202609052009
place $U/OS itstools 202609051237 OS
place $A/RD-2026-09-05 itstools 202609051237 RD
place $U/warmup-2026-09-06 itstools 202609060003
place $U/warmup-2026-09-07 itstools 2026090711
place $A/2026-09-07 itstools 202609071104
place $A/202609080313 itstools 202609080313
place $U/20260908-its itstools 2026090817
place $A/20260908-hsc hsc 20260908
mkdir -p "$R/_shared"
for x in oracle-2026 oracles-merged-2026 pages total support.csv pages-build.log pages-build-2026-09-07.log pages-build-202609080313.log pages-build-20260908-hsc.log pages-build-20260908-its.log; do [ -e "$U/$x" ] && mv "$U/$x" "$R/_shared/"; done
[ -e "$A/old" ] && mv "$A/old" "$R/_shared/old"
rmdir "$U" "$A" 2>/dev/null && echo "removed empty MCC26run and MCC26archive" || { echo "left over:"; ls "$U" "$A" 2>/dev/null; }
