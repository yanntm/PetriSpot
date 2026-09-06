#!/usr/bin/env python3
"""Does low coverage in a total examination follow the size of the net?

Reads the rows of one examination from a `total-runs.csv` (totallogs2csv.py),
pulls the place, transition and arc counts from the `Parsed PT model` line of
each log, and prints the runs grouped by transition count decile: how many
reached the wall, their median completion and median atoms answered, the
walker's share; then the wall runs with the fewest atoms answered.

    qla_coverage.py total-runs.csv [--exam QuasiLivenessAll] [--wall 1750]
"""

import argparse
import csv
import re
import statistics

RE_PARSED = re.compile(r"Parsed PT model containing (\d+) places and (\d+) transitions and (\d+) arcs")


def sizes(log):
    with open(log, errors="replace") as f:
        for line in f:
            m = RE_PARSED.search(line)
            if m:
                return int(m.group(1)), int(m.group(2)), int(m.group(3))
    return None, None, None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("csv")
    ap.add_argument("--exam", default="QuasiLivenessAll")
    ap.add_argument("--wall", type=float, default=1750, help="seconds above which a run is at the wall")
    args = ap.parse_args()
    rows = []
    with open(args.csv, newline="") as f:
        for r in csv.DictReader(f):
            if r["Examination"] != args.exam or not int(r["atoms"] or 0):
                continue
            p, t, a = sizes(r["log"])
            if t is None:
                continue
            rows.append({"model": r["Model"], "P": p, "T": t, "A": a, "atoms": int(r["atoms"]), "answered": int(r["answered"]),
                         "completion": float(r["completion"]), "s": int(r["duration(ms)"] or 0) / 1000,
                         "walk calls": int(r["walk calls"] or 0), "walk s": int(r["walk ms"] or 0) / 1000,
                         "wall": int(r["duration(ms)"] or 0) / 1000 > args.wall})
    rows.sort(key=lambda r: r["T"])
    n = len(rows)
    print(f"{args.exam}: {n} runs, {sum(r['wall'] for r in rows)} at the wall\n")
    print(f"{'transitions':>22s} {'runs':>5s} {'wall':>5s} {'median completion':>18s} {'median answered':>16s} {'median arcs/place':>18s} {'walker s (median)':>18s}")
    for d in range(10):
        chunk = rows[d * n // 10:(d + 1) * n // 10]
        if not chunk:
            continue
        wall = [r for r in chunk if r["wall"]]
        med = lambda xs: statistics.median(xs) if xs else float("nan")
        print(f"{chunk[0]['T']:>10d} .. {chunk[-1]['T']:<9d} {len(chunk):5d} {len(wall):5d} "
              f"{med([r['completion'] for r in chunk]):18.3f} {med([r['answered'] for r in chunk]):16.0f} "
              f"{med([r['A'] / max(r['P'], 1) for r in chunk]):18.1f} {med([r['walk s'] for r in chunk]):18.0f}")
    print("\nRuns at the wall, fewest atoms answered first (answered / atoms, transitions, arcs per place, walker calls and seconds):")
    for r in sorted((r for r in rows if r["wall"]), key=lambda r: r["answered"])[:40]:
        print(f"  {r['model']:38s} {r['answered']:7d} / {r['atoms']:<7d} T={r['T']:<7d} A/P={r['A'] / max(r['P'], 1):8.1f} walk {r['walk calls']} / {r['walk s']:.0f} s")
    wall = [r for r in rows if r["wall"]]
    if wall:
        big = [r for r in wall if r["T"] >= 20000]
        print(f"\nAt the wall: {len(wall)} runs; {len(big)} have at least 20000 transitions, "
              f"median completion {statistics.median([r['completion'] for r in big]) if big else float('nan'):.3f} against "
              f"{statistics.median([r['completion'] for r in wall if r['T'] < 20000]) if len(big) < len(wall) else float('nan'):.3f} for the others.")


if __name__ == "__main__":
    main()
