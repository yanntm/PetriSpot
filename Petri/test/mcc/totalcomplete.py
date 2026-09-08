#!/usr/bin/env python3
"""Complete the vector oracles where a one bit answer determines every object.

Two of the contest's ordinary examinations pin a whole vector:

* `QuasiLiveness` TRUE says *every* transition is quasi-live, so the QLA
  vector is all `T`.
* `StableMarking` FALSE says *no* place is stable, so the SMA vector is
  all `F`.

The other directions determine nothing (a FALSE QuasiLiveness says only that
some transition is not quasi-live), and `OneSafe` TRUE bounds every place by
one without saying which are ever marked, so it completes nothing either.

The source of those answers is the curated oracle of the corpus itself, the
`<model>-QL.out` and `<model>-SM.out` files beside the vectors: they are the
consensus to trust, ahead of `raw-result-analysis.csv`. A value already in a
vector is never overwritten: one that disagrees is reported and left alone,
since it means the oracle and our runs contradict each other and neither
should be quietly dropped.

    totalcomplete.py ORACLE_DIR [-o OUT_DIR]
"""
import argparse
import collections
import os
import re
from typing import Dict, List, Tuple

RE_VERDICT = re.compile(r"^FORMULA (\S+) (TRUE|FALSE)\b")

FILLS = {("QuasiLiveness", "TRUE"): ("QLA", "T"), ("StableMarking", "FALSE"): ("SMA", "F")}


def read_consensus(oracles: str) -> Dict[Tuple[str, str], str]:
    """(model, examination) -> TRUE/FALSE from the corpus's own oracle files."""
    out: Dict[Tuple[str, str], str] = {}
    for exam, suffix in (("QuasiLiveness", "QL"), ("StableMarking", "SM")):
        for name in sorted(os.listdir(oracles)):
            if not name.endswith(f"-{suffix}.out"):
                continue
            model = name[: -len(f"-{suffix}.out")]
            with open(os.path.join(oracles, name), errors="replace") as f:
                for line in f:
                    m = RE_VERDICT.match(line)
                    if m and m.group(1) == exam:
                        out[(model, exam)] = m.group(2)
                        break
    return out


def read_vector(path: str) -> Tuple[str, str, List[str]]:
    with open(path, errors="replace") as f:
        lines = [l.rstrip("\n") for l in f]
    return lines[0], lines[1].strip(), [c for line in lines[2:] for c in line.strip()]


def write_vector(path: str, header: str, keyword: str, values: List[str]) -> None:
    with open(path, "w") as f:
        f.write(f"{header}\n{keyword}\n")
        chars = "".join(values)
        for i in range(0, len(chars), 80):
            f.write(chars[i:i + 80] + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("oracles", help="directory of the oracle files: the <QLA|SMA> vectors and the <QL|SM> verdicts")
    ap.add_argument("-o", "--out", help="write here instead of in place")
    args = ap.parse_args()
    outdir = args.out or args.oracles
    if args.out:
        os.makedirs(outdir, exist_ok=True)

    consensus = read_consensus(args.oracles)
    stats: Dict[str, collections.Counter] = collections.defaultdict(collections.Counter)
    disagreements = 0
    for (model, exam), verdict in sorted(consensus.items()):
        fill = FILLS.get((exam, verdict))
        if fill is None:
            continue
        suffix, token = fill
        path = os.path.join(args.oracles, f"{model}-{suffix}.out")
        if not os.path.exists(path):
            stats[suffix]["no vector"] += 1
            continue
        header, keyword, values = read_vector(path)
        clash = [i for i, v in enumerate(values) if v != "?" and v != token]
        if clash:
            disagreements += 1
            stats[suffix]["DISAGREES with the consensus"] += 1
            print(f"DISAGREEMENT {model} {suffix}: the consensus says {exam} {verdict} "
                  f"(every object {token}) but the vector holds {values[clash[0]]} at {clash[0]}"
                  f" and {len(clash) - 1} more" if len(clash) > 1 else "")
            continue
        filled = sum(1 for v in values if v == "?")
        if filled:
            stats[suffix]["completed"] += 1
            stats[suffix]["objects filled"] += filled
            write_vector(os.path.join(outdir, f"{model}-{suffix}.out"), header, keyword,
                         [token] * len(values))
        else:
            stats[suffix]["already complete"] += 1
    for suffix in ("QLA", "SMA"):
        if stats[suffix]:
            print(f"{suffix}: " + ", ".join(f"{k} {v}" for k, v in sorted(stats[suffix].items())))
    return 1 if disagreements else 0


if __name__ == "__main__":
    raise SystemExit(main())
