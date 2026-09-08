#!/usr/bin/env python3
"""Merge the vector oracles of several runs into one, and report disagreement.

A total examination answers one value per object: `T`/`F` per transition for
QuasiLivenessAll, per place for StableMarkingAll, an integer per place for
UpperBoundsAll, and `?` wherever a run said nothing. Several runs of the same
model answer different subsets, so merging them covers more than any one run.

The merge is only sound if the runs agree, so that is what this reports
first: for every object where two runs give *different* known values, the
value is dropped to `?` and the disagreement is listed. A disagreement means
one of the runs is wrong, and no oracle should be published over it.

    mergeoracles.py DIR...                        # report only
    mergeoracles.py -o OUT DIR...                 # write the merged set
    mergeoracles.py -o OUT --base CONTEST DIR...  # and say what it adds to CONTEST

Each DIR holds `<model>-<QLA|SMA|UBA>.out` files (as `log2oracle.py -o` or
`totallogs2csv.py --oracles` write them).
"""
import argparse
import collections
import os
import sys
from typing import Dict, List, Optional, Tuple

SUFFIXES = ("QLA", "SMA", "UBA")


def read_oracle(path: str) -> Optional[Tuple[str, str, str, List[str]]]:
    """(model, examination, keyword, values) of one oracle file, or None."""
    with open(path, errors="replace") as f:
        lines = [l.rstrip("\n") for l in f]
    if len(lines) < 2:
        return None
    head = lines[0].split()
    if len(head) < 2:
        return None
    model, exam = head[0], head[1]
    keyword = lines[1].strip()
    body = lines[2:]
    if keyword == "BOUND":
        values = [t for line in body for t in line.split()]
    else:
        values = [c for line in body for c in line.strip()]
    return model, exam, keyword, values


def merge(values: List[List[str]]) -> Tuple[List[str], List[Tuple[int, List[str]]]]:
    """One vector from several, plus the positions where they disagree."""
    width = max(len(v) for v in values)
    out: List[str] = []
    clashes: List[Tuple[int, List[str]]] = []
    for i in range(width):
        known = {v[i] for v in values if i < len(v) and v[i] != "?"}
        if not known:
            out.append("?")
        elif len(known) == 1:
            out.append(known.pop())
        else:
            out.append("?")
            clashes.append((i, sorted(known)))
    return out, clashes


def write_oracle(path: str, model: str, exam: str, keyword: str, values: List[str]) -> None:
    """The merged vector, in the format the contest oracles use."""
    with open(path, "w") as f:
        f.write(f"{model} {exam}\n{keyword}\n")
        if keyword == "BOUND":
            line = ""
            for t in values:
                if len(line) + len(t) + 1 > 80:
                    f.write(line + "\n")
                    line = t
                else:
                    line = t if not line else line + " " + t
            f.write(line + "\n")
        else:
            chars = "".join(values)
            for i in range(0, len(chars), 80):
                f.write(chars[i:i + 80] + "\n")


def known(values: List[str]) -> int:
    return sum(1 for v in values if v != "?")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dirs", nargs="+", help="directories of <model>-<QLA|SMA|UBA>.out files")
    ap.add_argument("-o", "--outdir", help="write the merged oracles here")
    ap.add_argument("--base", help="a directory of oracles to compare the merge against")
    args = ap.parse_args()

    # gather: (model, suffix) -> [(dir, model, exam, keyword, values)]
    seen: Dict[Tuple[str, str], List[Tuple[str, str, str, str, List[str]]]] = collections.defaultdict(list)
    for d in args.dirs:
        if not os.path.isdir(d):
            print(f"not a directory, skipped: {d}", file=sys.stderr)
            continue
        for name in sorted(os.listdir(d)):
            suffix = name[-7:-4]
            if not name.endswith(".out") or suffix not in SUFFIXES:
                continue
            got = read_oracle(os.path.join(d, name))
            if got is None:
                print(f"unreadable oracle, skipped: {os.path.join(d, name)}", file=sys.stderr)
                continue
            model, exam, keyword, values = got
            seen[(model, suffix)].append((d, model, exam, keyword, values))

    if args.outdir:
        os.makedirs(args.outdir, exist_ok=True)

    stats = collections.defaultdict(lambda: [0, 0, 0, 0])  # suffix -> files, objects, known, clashes
    clashing_models: Dict[str, List[str]] = collections.defaultdict(list)
    for (model, suffix), entries in sorted(seen.items()):
        _, _, exam, keyword, _ = entries[0]
        values, clashes = merge([e[4] for e in entries])
        s = stats[suffix]
        s[0] += 1
        s[1] += len(values)
        s[2] += known(values)
        s[3] += len(clashes)
        if clashes:
            clashing_models[suffix].append(model)
            for i, vals in clashes[:4]:
                sources = ", ".join(f"{os.path.basename(e[0])}={e[4][i] if i < len(e[4]) else '?'}" for e in entries)
                print(f"DISAGREEMENT {model} {suffix} object {i}: {vals} ({sources})")
        if args.outdir:
            write_oracle(os.path.join(args.outdir, f"{model}-{suffix}.out"), model, exam, keyword, values)

    print()
    for suffix in SUFFIXES:
        files, objects, k, clashes = stats[suffix]
        if not files:
            continue
        base_known = 0
        if args.base:
            for (model, sfx) in seen:
                if sfx != suffix:
                    continue
                p = os.path.join(args.base, f"{model}-{suffix}.out")
                got = read_oracle(p) if os.path.exists(p) else None
                if got:
                    base_known += known(got[3])
        share = 100.0 * k / objects if objects else 0.0
        line = f"{suffix}: {files} models, {objects} objects, {k} known ({share:.1f}%), {clashes} disagreements"
        if args.base:
            line += f"; the base knows {base_known}"
        print(line)
        if clashing_models[suffix]:
            print(f"    models with a disagreement: {' '.join(sorted(clashing_models[suffix])[:8])}")
    return 1 if any(s[3] for s in stats.values()) else 0


if __name__ == "__main__":
    raise SystemExit(main())
