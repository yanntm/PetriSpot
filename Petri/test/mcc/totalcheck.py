#!/usr/bin/env python3
"""Cross-check the vector oracles of the total examinations with the contest.

A model is QuasiLive iff every transition is quasi-live, and has a StableMarking
iff some place is stable. So a QLA vector holding an `F` says QuasiLiveness is
FALSE, an all-`T` vector says TRUE; an SMA vector holding a `T` says
StableMarking is TRUE, an all-`F` vector says FALSE. The contest consensus for
those two examinations (`estimated result` of `raw-result-analysis.csv`) is
the only independent check the total vectors have; this compares the two and
reports every contradiction, every confirmation, and the vectors the
consensus cannot decide.

    totalcheck.py raw-result-analysis.csv ORACLE_DIR [-o totalcheck.csv]

ORACLE_DIR holds `<model>-QLA.out` / `-SMA.out` files as written by
`log2oracle.py` or `totallogs2csv.py --oracles`.
"""

import argparse
import collections
import csv
import os
import sys

import toolsupport as ts


def read_vector(path):
    """The verdict characters of a QLA/SMA vector (T, F, ?)."""
    with open(path) as f:
        lines = f.read().split("\n")
    return "".join(lines[2:]).replace(" ", "")


def implied(vec, exam):
    """The global verdict a vector implies, or None when the vector cannot decide."""
    if exam == "QLA":
        if "F" in vec:
            return "FALSE"
        if "?" not in vec:
            return "TRUE"
    else:
        if "T" in vec:
            return "TRUE"
        if "?" not in vec:
            return "FALSE"
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("raw", help="raw-result-analysis.csv")
    ap.add_argument("oracles", help="directory of <model>-QLA.out / -SMA.out vectors")
    ap.add_argument("-o", "--out", help="write one row per vector as CSV")
    args = ap.parse_args()

    consensus = {}
    with open(args.raw, newline="") as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            if row[0] in ts.VIRTUAL or row[2] not in ("QuasiLiveness", "StableMarking"):
                continue
            est = row[15].strip()
            consensus[(row[1], "QLA" if row[2] == "QuasiLiveness" else "SMA")] = ts.normalize(est) if est not in ts.NO_ANSWER else "?"

    census = collections.Counter()
    rows = []
    for name in sorted(os.listdir(args.oracles)):
        for exam in ("QLA", "SMA"):
            if name.endswith(f"-{exam}.out"):
                model = name[: -len(f"-{exam}.out")]
                vec = read_vector(os.path.join(args.oracles, name))
                ours = implied(vec, exam)
                theirs = consensus.get((model, exam), "?")
                if ours is None:
                    status = "undecided by the vector"
                elif theirs == "?":
                    status = "consensus unknown"
                elif ours == theirs:
                    status = "confirmed"
                else:
                    status = "CONTRADICTION"
                    print(f"CONTRADICTION {model} {exam}: vector implies {ours}, consensus {theirs} "
                          f"({vec.count('T')} T, {vec.count('F')} F, {vec.count('?')} ? of {len(vec)})")
                census[(exam, status)] += 1
                rows.append({"Model": model, "Examination": exam, "atoms": len(vec), "T": vec.count("T"),
                             "F": vec.count("F"), "?": vec.count("?"), "implied": ours or "", "consensus": theirs, "status": status})
    for (exam, status), n in sorted(census.items()):
        print(f"{exam} {status:24s} {n}", file=sys.stderr)
    if args.out:
        with open(args.out, "w", newline="") as f:
            w = csv.DictWriter(f, ["Model", "Examination", "atoms", "T", "F", "?", "implied", "consensus", "status"])
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    main()
