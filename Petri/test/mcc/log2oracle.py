#!/usr/bin/env python3
"""Extract an oracle file from the harness log of a total examination run.

The log holds the tool's answers, `QLIVE t12 TRUE ...`, `BOUND p3 7 ...`; the
oracle is the same vector in the pnmcc-models-2026 format, the header line,
the keyword, then `T`/`F`/`?` per object wrapped at 80 columns (one token per
place for the bounds), a `?` wherever the run said nothing.

    log2oracle.py OAR.123.stdout                 # the oracle on stdout
    log2oracle.py -o DIR QLA/*.stdout            # one <model>-QLA.out per log in DIR
"""

import argparse
import os
import sys

import totallogs2csv as total


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("logs", nargs="+", help="harness logs of total examination runs")
    ap.add_argument("-o", "--outdir", help="write <model>-<QLA|SMA|UBA>.out files here instead of stdout")
    args = ap.parse_args()

    if args.outdir:
        os.makedirs(args.outdir, exist_ok=True)
    for path in args.logs:
        run, atoms, verdict = total.parse(path)
        if run["Examination"] not in total.SUFFIX or not atoms:
            print(f"{path}: not a total examination log", file=sys.stderr)
            continue
        if args.outdir:
            total.write_oracle(args.outdir, run, atoms, verdict)
        else:
            total.write_oracle_to(sys.stdout, run, atoms, verdict)


if __name__ == "__main__":
    main()
