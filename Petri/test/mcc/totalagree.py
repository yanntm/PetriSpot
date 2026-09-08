#!/usr/bin/env python3
"""Do our own vector oracles agree with our own single-verdict answers?

A total examination answers one value per object; three ordinary examinations
answer the same question in one bit, and the two must agree:

* `QuasiLiveness` is TRUE exactly when every transition is quasi-live, so a
  TRUE forbids an `F` in the QLA vector and a FALSE demands one.
* `StableMarking` is TRUE exactly when some place is stable, so a TRUE demands
  a `T` in the SMA vector and a FALSE forbids one.
* `OneSafe` is TRUE exactly when no place ever holds more than one token, so a
  TRUE forbids a bound above 1 in the UBA vector and a FALSE demands one.

A vector with unknowns can only contradict in one direction, which is why each
case is counted apart: confirmed, contradicted, or undecided by the vector.
This is an internal check — the same tool answering the same question two
ways — and a contradiction is a bug in one of the two.

    totalagree.py ORACLE_DIR --ql QL_DIR --sm SM_DIR --os OS_DIR
"""
import argparse
import collections
import os
import re
import sys
from typing import Dict, List, Optional, Tuple

RE_MODEL = re.compile(r"runatest[_a-z]*\.sh ([A-Za-z0-9_.-]+)")
RE_FORMULA = re.compile(r"^FORMULA (\S+) (\S+)")


def verdicts_of(directory: str, prop: str) -> Dict[str, str]:
    """model -> TRUE/FALSE for the named property, from harness logs."""
    out: Dict[str, str] = {}
    if not directory or not os.path.isdir(directory):
        return out
    for name in sorted(os.listdir(directory)):
        if not (name.endswith("out") or name.endswith(".log")):
            continue
        model = None
        value = None
        with open(os.path.join(directory, name), errors="replace") as f:
            for line in f:
                if model is None:
                    m = RE_MODEL.search(line)
                    if m:
                        model = m.group(1)
                m = RE_FORMULA.match(line)
                if m and m.group(1) == prop:
                    value = m.group(2)
        if model and value in ("TRUE", "FALSE"):
            out[model] = value
    return out


def read_vector(path: str) -> Optional[Tuple[str, List[str]]]:
    """(keyword, values) of a vector oracle."""
    if not os.path.exists(path):
        return None
    with open(path, errors="replace") as f:
        lines = [l.rstrip("\n") for l in f]
    if len(lines) < 2:
        return None
    keyword = lines[1].strip()
    body = lines[2:]
    if keyword == "BOUND":
        return keyword, [t for line in body for t in line.split()]
    return keyword, [c for line in body for c in line.strip()]


def bound_above_one(token: str) -> Optional[bool]:
    """Is this bound above one? None when unknown."""
    if token == "?":
        return None
    if token in ("inf", "INF", "oo"):
        return True
    try:
        return int(token) > 1
    except ValueError:
        return None


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("oracles", help="directory of <model>-<QLA|SMA|UBA>.out vectors")
    ap.add_argument("--ql", help="directory of QuasiLiveness logs")
    ap.add_argument("--sm", help="directory of StableMarking logs")
    ap.add_argument("--os", dest="onesafe", help="directory of OneSafe logs")
    args = ap.parse_args()

    cases = [("QuasiLiveness", args.ql, "QLA"), ("StableMarking", args.sm, "SMA"),
             ("OneSafe", args.onesafe, "UBA")]
    bad = 0
    for prop, directory, suffix in cases:
        verdicts = verdicts_of(directory, prop)
        if not verdicts:
            print(f"{prop}: no verdicts read")
            continue
        c: collections.Counter = collections.Counter()
        for model, verdict in sorted(verdicts.items()):
            got = read_vector(os.path.join(args.oracles, f"{model}-{suffix}.out"))
            if got is None:
                c["no vector"] += 1
                continue
            _, values = got
            if suffix == "QLA":       # TRUE: every transition quasi-live
                witness = any(v == "F" for v in values)      # forbids TRUE
                complete = all(v != "?" for v in values)
            elif suffix == "SMA":     # TRUE: some place stable
                witness = any(v == "T" for v in values)      # demands TRUE
                complete = all(v != "?" for v in values)
            else:                     # OneSafe TRUE: no bound above one
                flags = [bound_above_one(v) for v in values]
                witness = any(f is True for f in flags)      # forbids TRUE
                complete = all(f is not None for f in flags)
            demands_true = witness if suffix == "SMA" else not witness
            if verdict == "TRUE":
                if (suffix == "SMA" and witness) or (suffix != "SMA" and not witness and complete):
                    c["confirmed"] += 1
                elif (suffix == "SMA" and not witness and complete) or (suffix != "SMA" and witness):
                    c["CONTRADICTED"] += 1
                    print(f"CONTRADICTION {model}: {prop} TRUE against its {suffix} vector")
                    bad += 1
                else:
                    c["undecided by the vector"] += 1
            else:  # FALSE
                if (suffix == "SMA" and not witness and complete) or (suffix != "SMA" and witness):
                    c["confirmed"] += 1
                elif (suffix == "SMA" and witness) or (suffix != "SMA" and not witness and complete):
                    c["CONTRADICTED"] += 1
                    print(f"CONTRADICTION {model}: {prop} FALSE against its {suffix} vector")
                    bad += 1
                else:
                    c["undecided by the vector"] += 1
            _ = demands_true
        print(f"{prop} against {suffix}: " + ", ".join(f"{k} {v}" for k, v in sorted(c.items())))
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
