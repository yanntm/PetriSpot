#!/usr/bin/env python3
"""The exception behind each "eclipse fatal error" of a campaign.

When the ITS-Tools application thread dies on an uncaught exception, the
Eclipse launcher tries to show a dialog, fails without a display, and prints
`An error has occurred. See the log file` on stdout; the exception itself goes
to `itstools/configuration/<epoch ms>.log` on the cluster. This joins the two:
each Eclipse log is matched to the run whose ITS-Tools start stamp is nearest
to the log's session time, among the runs the collector marked `eclipse_fatal`.

    eclipselogs.py ECLIPSE_LOG_DIR runs.csv [runs.csv ...] [-o eclipse.csv]

Prints one line per Eclipse log (session time, exception, matched run) and a
census of the exceptions.
"""

import argparse
import collections
import csv
import datetime
import os
import re
import sys

RE_SESSION = re.compile(r"^!SESSION (\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)")
RE_EXC = re.compile(r"^([a-zA-Z_.$]+(?:Exception|Error))(?::\s*(.*))?$")


def read_eclipse_log(path):
    """(session time, exception, message) of one configuration log."""
    session = exc = msg = ""
    with open(path, errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            m = RE_SESSION.match(line)
            if m and not session:
                session = m.group(1)
            m = RE_EXC.match(line)
            if m and not exc:
                exc, msg = m.group(1), (m.group(2) or "")[:100]
    return session, exc, msg


def read_runs(paths):
    """The eclipse_fatal runs: (start time, model, examination, log)."""
    runs = []
    for p in paths:
        with open(p, newline="") as f:
            for r in csv.DictReader(f):
                if r["failure"] == "eclipse_fatal" and r["start"]:
                    runs.append((r["start"], r["Model"], r["Examination"], r["log"]))
    return runs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("logdir", help="directory of Eclipse configuration/*.log files")
    ap.add_argument("runs", nargs="+", help="runs.csv files of the campaigns")
    ap.add_argument("-o", "--out", help="write the join as CSV")
    args = ap.parse_args()

    runs = read_runs(args.runs)
    stamps = [(datetime.datetime.strptime(s, "%Y-%m-%d %H:%M:%S"), m, e, l) for s, m, e, l in runs]
    rows = []
    census = collections.Counter()
    for name in sorted(os.listdir(args.logdir)):
        if not name.endswith(".log"):
            continue
        session, exc, msg = read_eclipse_log(os.path.join(args.logdir, name))
        match = ("", "", "")
        if session and stamps:
            t = datetime.datetime.strptime(session, "%Y-%m-%d %H:%M:%S")
            # the run starts before the Eclipse session; take the nearest one within a day
            cands = [(t - s, m, e, l) for s, m, e, l in stamps if datetime.timedelta(0) <= t - s < datetime.timedelta(days=1)]
            if cands:
                d, m, e, l = min(cands)
                match = (m, e, f"{d.total_seconds():.0f}s after start")
        census[exc or "(no exception line)"] += 1
        rows.append({"eclipse log": name, "session": session, "exception": exc, "message": msg,
                     "Model": match[0], "Examination": match[1], "offset": match[2]})
        print(f"{session}  {exc:45s} {match[0]:40s} {match[1]:24s} {match[2]}")

    print("\nexceptions:", file=sys.stderr)
    for k, v in census.most_common():
        print(f"  {v:4d}  {k}", file=sys.stderr)
    print(f"{len(rows)} Eclipse logs, {len(runs)} eclipse_fatal runs in the CSVs", file=sys.stderr)
    if args.out:
        with open(args.out, "w", newline="") as f:
            w = csv.DictWriter(f, ["eclipse log", "session", "exception", "message", "Model", "Examination", "offset"])
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    main()
