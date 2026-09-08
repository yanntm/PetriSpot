#! /bin/bash
# What the Louvain decomposition is handed and what it costs, on a ladder of instances.
#
#   PROD=/data/ythierry/itstools-local bash Petri/test/louvain-bench.sh [model...]
#
# Needs a local product built with the -louvainBench flag (Application) and GraphBuilder
# DEBUG >= 2: the flag parses the model and its properties, runs the decomposition alone --
# no engine, no solver, so nothing answers the properties first and the measurement is
# deterministic -- and the DEBUG level keeps the graph, binary, weights and tree of each
# call and prints their paths, which is what this reads.
#
# Per model: the variables and constraints handed to the graph, the edges those induce, the
# graph written, and the time convert and louvain need on it. Logs in Petri/test/logs/louvain-bench/.
set -u
PROD=${PROD:-/data/ythierry/itstools-local}
IN=${IN:-/data/ythierry/MCC26deploy/MCC-drivers/INPUTS}
EXAM=${EXAM:-ReachabilityCardinality}
TIMEOUT=${TIMEOUT:-60}
OUT=${OUT:-Petri/test/logs/louvain-bench}
BIN=$(ls -d "$PROD"/plugins/fr.lip6.move.gal.louvain.binaries_*/bin)
mkdir -p "$OUT"
models=${*:-"Philosophers-PT-000005 Philosophers-PT-000010 Philosophers-PT-000020 Philosophers-PT-000050 Philosophers-PT-000100 SharedMemory-PT-000005 SharedMemory-PT-000010 SharedMemory-PT-000050 AirplaneLD-PT-0010"}
printf "%-26s %6s %6s %7s %10s %10s %8s %8s %6s %7s\n" model vars cstrs maxcstr edges bytes convert louvain comms decomp_ms
for m in $models ; do
	[ -d "$IN/$m" ] || tar xzf "$IN/$m.tgz" -C "$IN" 2>/dev/null
	log="$OUT/$m.log"
	rm -f /tmp/graph*
	timeout $((TIMEOUT * 3)) "$PROD"/its-tools -pnfolder "$IN/$m" -examination "$EXAM" -louvainBench -order META -timeout "$TIMEOUT" > "$log" 2>&1
	vars=$(sed -n 's/^Louvain bench : \([0-9]*\) variables,.*/\1/p' "$log" | head -1)
	cstrs=$(sed -n 's/^Louvain bench : \([0-9]*\) constraints inducing.*/\1/p' "$log" | head -1)
	maxc=$(sed -n 's/^Louvain bench : .* support [0-9]* median \([0-9]*\) max.*/\1/p' "$log" | head -1)
	decomp=$(sed -n 's/^Louvain bench : decomposition took \([0-9]*\) ms/\1/p' "$log" | tail -1)
	edges=0 ; bytes=0 ; conv=0 ; louv=0 ; comms=0
	for txt in /tmp/graph*.txt ; do
		[ -f "$txt" ] || continue
		edges=$((edges + $(wc -l < "$txt")))
		bytes=$((bytes + $(stat -c%s "$txt")))
		c=$( { /usr/bin/time -f %e "$BIN"/convert-linux64 -i "$txt" -o /tmp/bench.bin -w /tmp/bench.w ; } 2>&1 | tail -1)
		l=$( { /usr/bin/time -f %e "$BIN"/louvain-linux64 /tmp/bench.bin -l -1 -v -w /tmp/bench.w -q 0 -e 0.001 > /tmp/bench.tree ; } 2>&1 | tail -1)
		conv=$(echo "$conv + $c" | bc) ; louv=$(echo "$louv + $l" | bc)
		comms=$(awk '{print $2}' /tmp/bench.tree | sort -u | wc -l)
	done
	printf "%-26s %6s %6s %7s %10d %10d %8s %8s %6d %7s\n" "$m" "${vars:--}" "${cstrs:--}" "${maxc:--}" \
	  "$edges" "$bytes" "$conv" "$louv" "$comms" "${decomp:--}"
	rm -f /tmp/bench.bin /tmp/bench.w /tmp/bench.tree
done
