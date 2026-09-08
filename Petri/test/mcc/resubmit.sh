#! /bin/bash
# resubmit.sh: submit the instances of an examination that have no usable log yet.
#
#   ssh cluster.lip6.fr 'cd ~/MCC26/MCC-drivers && TIMEOUT=1800 WALLTIME=0:45:00 \
#      HOSTS=tall% CORES=4 RUNATEST=./runatest_cluster.sh ./resubmit.sh CTLC'
#
# Runs on the cluster head, from the harness root. An instance is redone when it
# has no log, when its log does not reach the end, or when the log carries one of
# the signatures in BAD: a native image closed world miss answers nothing in a
# second, so such a log is worth nothing and costs nothing to replace. Everything
# else is left alone, which is what makes an interrupted campaign cheap to finish.
#
# A run is complete only when it prints the teamcity suite close. "Running test"
# and "Test :" both appear at the start, so neither tells a finished run from one
# killed at the wall or by oardel -- of the 177 logs an interrupted CTLC left
# behind, all 177 carried "Running test" and only 67 had actually finished.
set -eu
EX=${1:?examination folder, e.g. CTLC}
BAD=${BAD:-MissingReflectionRegistrationError}
TREE=$(cd "$(dirname "$0")" && pwd)
cd "$TREE"
tmp=$(mktemp -d "$TREE/.resubmit.XXXXXX")
trap 'rm -rf "$tmp"' EXIT

ls oracle/*-"$EX".out | sed "s|^oracle/||; s|-$EX\.out$||" | sort > "$tmp/all"
: > "$tmp/done"
if [ -d "$EX" ]; then
	for f in "$EX"/*.stdout ; do
		[ -e "$f" ] || continue
		grep -qE "$BAD" "$f" && continue
		grep -q "testSuiteFinished" "$f" || continue
		sed -n "s/.*testSuiteFinished name='oracle\.\(.*\)-$EX'.*/\1/p" "$f" | head -1 >> "$tmp/done"
	done
fi
sort -u "$tmp/done" -o "$tmp/done"
comm -23 "$tmp/all" "$tmp/done" | sed "s|^|oracle/|; s|\$|-$EX.out|" > "$tmp/todo"

n=$(wc -l < "$tmp/todo")
echo "$EX: $(wc -l < "$tmp/all") instances, $(wc -l < "$tmp/done") already good, $n to submit"
[ "$n" -gt 0 ] || exit 0
[ "${DRYRUN:-0}" = 1 ] && { head -5 "$tmp/todo"; echo "(dry run)"; exit 0; }
./run_oar.sh "$(tr '\n' ' ' < "$tmp/todo")"
