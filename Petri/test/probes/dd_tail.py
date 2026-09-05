#!/usr/bin/env python3
"""How long does a run spend after its last word from the Java side?

ITS-Tools timestamps its own `[INFO]` lines; the decision diagram binaries do
not. So the gap between the last timestamped line and the end of the suite is
the tail in which the portfolio has nothing left to say and only the (single
threaded) diagram engine is still working. That tail is the window in which
other cores are available.

    dd_tail.py runs.csv <logdir>...
"""

import csv
import os
import re
import sys
from datetime import datetime

RE_STAMP = re.compile(r"^\[(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)\]")
FMT = "%Y-%m-%d %H:%M:%S"


def span(path):
    first = last = None
    for line in open(path, errors="replace"):
        m = RE_STAMP.match(line)
        if m:
            if first is None:
                first = m.group(1)
            last = m.group(1)
    if first is None or last is None:
        return None
    return (datetime.strptime(last, FMT) - datetime.strptime(first, FMT)).total_seconds()


def main(runs_csv, dirs):
    rows = {}
    with open(runs_csv) as f:
        for r in csv.DictReader(f):
            rows[r["log"]] = r
    by_exam = {}
    for d in dirs:
        for name in sorted(os.listdir(d)):
            if not name.endswith("out"):
                continue
            path = os.path.join(d, name)
            r = rows.get(path)
            if r is None or not r["duration(ms)"]:
                continue
            dur = int(r["duration(ms)"]) / 1000
            if dur < 1700:      # only the runs that used the whole budget
                continue
            s = span(path)
            if s is None:
                continue
            tail = max(0.0, dur - s)
            by_exam.setdefault(r["Examination"], []).append((tail, dur, float(r["walk ms"]) / 1000))

    print(f"{'examination':<22} {'runs':>5} {'mean tail':>10} {'share':>7} {'walk s':>8}")
    for exam, v in sorted(by_exam.items()):
        tails = [t for t, _, _ in v]
        share = sum(tails) / sum(d for _, d, _ in v)
        walk = sum(w for _, _, w in v) / len(v)
        print(f"{exam:<22} {len(v):>5} {sum(tails)/len(tails):>9.0f}s {share:>6.0%} {walk:>7.0f}s")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
