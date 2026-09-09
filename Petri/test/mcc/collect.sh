#!/usr/bin/env bash
# collect.sh: bring a campaign's logs down from the cluster and turn them into tables.
#
#   bash Petri/test/mcc/collect.sh <tool>/<build> EXAM... [--pages]
#   BASELINE=<csv folder> bash Petri/test/mcc/collect.sh itstools/202609080313 LTLC LTLF
#
# Rsyncs each examination folder of ~/MCC26/MCC-drivers on the cluster into the one log
# root, /data/ythierry/MCC26logs/<tool>/<build>/ (its README.md: the tool, then the
# product build id as the campaign label; never with --delete: a partial campaign is
# collected as it drains, and rerunning adds the new logs), then runs the collectors over
# every examination folder present there: mcclogs2csv.py for the classic examinations,
# totallogs2csv.py for QLA SMA UBA (the total examinations), and report.py, against
# $BASELINE when set. With --pages the MCC-analysis pages are rebuilt afterwards (about
# four minutes; the set whose glob names <tool>/<build>/* in campaign/example.json picks
# the new logs up; a new campaign needs a new set there first).
set -eu
HOST=${CLUSTER_HOST:-cluster.lip6.fr}
TREE=${CLUSTER_TREE:-MCC26/MCC-drivers}
RUNS=${RUNS:-/data/ythierry/MCC26logs}
MCC=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ANALYSIS=${ANALYSIS:-$HOME/git/MCC-analysis/campaign}
date=${1:?campaign, as <tool>/<build>, e.g. itstools/202609080313}; shift
pages=0; exams=""
for a in "$@"; do case $a in --pages) pages=1;; *) exams="$exams $a";; esac; done
[ -n "$exams" ] || { echo "no examination given"; exit 1; }
dest=$RUNS/$date; mkdir -p "$dest"
# report.py runs from $dest, so a relative BASELINE must be resolved here
if [ -n "${BASELINE:-}" ]; then BASELINE=$(cd "$BASELINE" && pwd); fi
for ex in $exams; do
  echo "== rsync $ex"
  rsync -rz --exclude='*.stderr' "$HOST:$TREE/$ex" "$dest/" || echo "rsync of $ex failed (folder absent on the cluster?)"
done
cd "$dest"
classic=""; total=""
for d in */; do d=${d%/}; [ "$d" = csv ] && continue; ls "$d"/*.stdout >/dev/null 2>&1 || continue
  case $d in QLA|SMA|UBA) total="$total $d";; *) classic="$classic $d";; esac; done
mkdir -p csv
[ -n "$classic" ] && { echo "== mcclogs2csv$classic"; python3 "$MCC/mcclogs2csv.py" $classic -o csv; }
[ -n "$total" ] && { echo "== totallogs2csv$total"; python3 "$MCC/totallogs2csv.py" $total -o csv --oracles csv/oracles; }
echo "== report"
if [ -n "${BASELINE:-}" ]; then python3 "$MCC/report.py" csv --baseline "$BASELINE"; else python3 "$MCC/report.py" csv; fi
echo "report: $dest/csv/REPORT.md"
if [ $pages = 1 ]; then
  echo "== pages (about four minutes)"
  python3 "$ANALYSIS/build.py" "$ANALYSIS/example.json" > "$RUNS/_shared/pages-build-${date//\//-}.log" 2>&1 && tail -2 "$RUNS/_shared/pages-build-${date//\//-}.log"
  echo "pages: $RUNS/_shared/pages (python3 $ANALYSIS/serve.py $RUNS/_shared/pages --port 8080 to browse)"
fi
