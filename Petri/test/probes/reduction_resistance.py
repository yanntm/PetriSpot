#!/usr/bin/env python3
"""Are the models that resist structural reductions the ones we fail on?

Every reduction pass reports

    Applied a total of N rules in M ms. Remains X /Y variables (removed Z) ...

This takes the first such line of each log as the yield of the first pass --
the fraction of places the reducer could remove before any property was solved
-- and crosses it with what the run ended up producing, read from a runs.csv
written by mcclogs2csv.py.

    reduction_resistance.py runs.csv <logdir>...
"""

import csv
import os
import re
import sys

RE_REMAINS = re.compile(r"Remains (\d+) */(\d+) variables \(removed (\d+)\)")


def first_yield(path):
    for line in open(path, errors="replace"):
        if line.startswith("Applied a total of"):
            m = RE_REMAINS.search(line)
            if m:
                total = int(m.group(2))
                return int(m.group(3)) / total if total else None
    return None


def main(runs_csv, dirs):
    outcome = {}
    with open(runs_csv) as f:
        for r in csv.DictReader(f):
            outcome[r["log"]] = r
    rows = []
    for d in dirs:
        for name in sorted(os.listdir(d)):
            if not name.endswith("out"):
                continue
            path = os.path.join(d, name)
            r = outcome.get(path)
            if r is None:
                continue
            y = first_yield(path)
            if y is None:
                continue
            rows.append((r["Examination"], r["Model"], y,
                         int(r["missed"]) + int(r["none"]), int(r["known"])))

    print(f"{len(rows)} logs reported a reduction pass\n")
    buckets = [(0.0, 0.01), (0.01, 0.1), (0.1, 0.3), (0.3, 0.6), (0.6, 1.01)]
    print(f"{'places removed by pass 1':<26} {'logs':>6} {'unanswered':>11} {'per log':>9}")
    for lo, hi in buckets:
        sel = [r for r in rows if lo <= r[2] < hi]
        if not sel:
            continue
        un = sum(r[3] for r in sel)
        print(f"{lo:>10.0%} - {hi:<12.0%} {len(sel):>6} {un:>11} {un/len(sel):>9.2f}")

    print()
    silent = [r for r in rows if r[3] > 0]
    speaking = [r for r in rows if r[3] == 0]
    for label, sel in (("runs that answered everything", speaking),
                       ("runs that left a value unanswered", silent)):
        if sel:
            mean = sum(r[2] for r in sel) / len(sel)
            zero = sum(1 for r in sel if r[2] < 0.01)
            print(f"{label:<36} {len(sel):>6} logs, mean yield {mean:>6.1%}, "
                  f"{100*zero/len(sel):>5.1f} % removed nothing")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
