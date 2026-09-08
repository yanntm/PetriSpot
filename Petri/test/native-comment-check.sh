#! /bin/bash
# One-shot: does the rebuilt image parse the COL models that died on OMCommentImpl?
# Runs each trouble family and AirplaneLD PT as a control, prints verdicts and any exception.
set -u
IMG=${IMG:-/data/ythierry/MCC26deploy/its-tools-native-test}
PROD=${PROD:-/data/ythierry/MCC26deploy/MCC-drivers/itstools/itstools}
IN=${IN:-/data/ythierry/MCC26deploy/MCC-drivers/INPUTS}
OUT=${OUT:-/home/ythierry/git/PetriSpot/Petri/test/logs/native-comment-check}
mkdir -p "$OUT"
for m in Sudoku-COL-AN01 LastZero-COL-N08 FileSystem-COL-N02I05B05 UtilityControlRoom-COL-Z2T4N02 VehicularWifi-COL-none AirplaneLD-PT-0010 ; do
	log="$OUT/$m.log"
	"$IMG" -Dfr.lip6.binaries.root="$PROD/plugins" -pnfolder "$IN/$m" -examination CTLCardinality \
	  -its -ltsmin -greatspnpath "$PROD/../greatspn/" -order META -manyOrder -smt -timeout 60 > "$log" 2>&1
	printf "%-32s %2d verdicts  %s\n" "$m" "$(grep -c '^FORMULA' "$log")" \
	  "$(grep -m1 -o 'MissingReflectionRegistrationError\|ClassNotFoundException\|NoSuchMethodException' "$log" || echo ok)"
done
