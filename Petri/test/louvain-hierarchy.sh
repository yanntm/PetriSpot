#! /bin/bash
# Does the decomposition survive, and what does the run answer?
#
#   PROD=/data/ythierry/MCC26deploy/products/itstools-local bash Petri/test/louvain-hierarchy.sh [model...]
#
# Per model and examination, through -louvainBench (the decomposition on its own, then the
# symbolic engine on what it produced): the GAL types the hierarchy holds -- one means it was
# fused back, libits reading a sum or a comparison within one component -- the time the
# decomposition took, and the verdicts the run answered. Logs in Petri/test/logs/louvain-hierarchy/.
set -u
PROD=${PROD:-/data/ythierry/MCC26deploy/products/itstools-local}
IN=${IN:-/data/ythierry/MCC26deploy/MCC-drivers/INPUTS}
LIMIT=${LIMIT:-/data/ythierry/MCC26deploy/MCC-drivers/limit_time.pl}
TIMEOUT=${TIMEOUT:-60}
EXAMS=${EXAMS:-"ReachabilityFireability ReachabilityCardinality"}
OUT=${OUT:-Petri/test/logs/louvain-hierarchy}
mkdir -p "$OUT"
OUT=$(cd "$OUT" && pwd)   # the run happens in the model folder
models=${*:-"AirplaneLD-PT-0010 Philosophers-PT-000020 SharedMemory-PT-000010 IOTPpurchase-PT-C05M04P03D02"}
printf "%-30s %-24s %14s %9s %8s\n" model examination hierarchy decomp_ms verdicts
for m in $models ; do
	[ -d "$IN/$m" ] || tar xzf "$IN/$m.tgz" -C "$IN" 2>/dev/null
	for e in $EXAMS ; do
		log="$OUT/$m-$e.log"
		(cd "$IN" && perl "$LIMIT" $((TIMEOUT * 3)) "$PROD"/its-tools -pnfolder "$m" -examination "$e" \
			-louvainBench -its -order META -timeout "$TIMEOUT" > "$log" 2>&1)
		types=$(sed -n 's/.*ms, \([0-9]*\) GAL types.*/\1/p' "$log" | tail -1)
		[ -z "$types" ] && grep -q "one flat GAL type" "$log" && types=1
		ms=$(sed -n 's/.*decomposition took \([0-9]*\) ms.*/\1/p' "$log" | tail -1)
		printf "%-30s %-24s %14s %9s %8s\n" "$m" "$e" "${types:--}" "${ms:--}" "$(grep -c '^FORMULA' "$log")"
	done
done
