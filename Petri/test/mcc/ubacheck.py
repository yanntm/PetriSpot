#!/usr/bin/env python3
"""Check the UpperBoundsAll vectors against the contest's UpperBounds consensus.

The sixteen UpperBounds formulas of a model bound a sum of places, often a
single place. For every single-place formula whose consensus value is known,
the UBA vector holds a token for that place: an integer must match, a `?` is
a miss, a different integer is a contradiction. Places are named by definition
index in the vector, so the PNML is read for the order of `<place id=...>`
and the UpperBounds.xml for the place of each formula; both come out of the
contest archives `<model>.tgz`.

    ubacheck.py INPUTS_DIR raw-result-analysis.csv ORACLE_DIR [-o ubacheck.csv]
"""

import argparse
import collections
import csv
import os
import re
import sys
import tarfile

import toolsupport as ts

RE_PLACE = re.compile(rb'<place\s[^>]*?id="([^"]+)"')
RE_PROP = re.compile(r"<property>(.*?)</property>", re.S)
RE_ID = re.compile(r"<id>(.*?)</id>")
RE_PLACEREF = re.compile(r"<place>(.*?)</place>")
RE_INDEX = re.compile(r"-(\d\d)$")


def read_archive(path, model):
    """The place ids in definition order, and {formula index: place id} of the single-place bounds."""
    places, single = [], {}
    with tarfile.open(path) as tar:
        for member in tar:
            base = os.path.basename(member.name)
            if base == "model.pnml":
                data = tar.extractfile(member).read()
                places = [m.decode() for m in RE_PLACE.findall(data)]
            elif base == "UpperBounds.xml":
                text = tar.extractfile(member).read().decode(errors="replace")
                for prop in RE_PROP.findall(text):
                    pid = RE_ID.search(prop)
                    refs = RE_PLACEREF.findall(prop)
                    m = RE_INDEX.search(pid.group(1).strip()) if pid else None
                    if m and len(refs) == 1:
                        single[int(m.group(1))] = refs[0].strip()
    return places, single


def read_vector(path):
    with open(path) as f:
        lines = f.read().split("\n")
    return " ".join(lines[2:]).split()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("inputs", help="directory of <model>.tgz archives")
    ap.add_argument("raw", help="raw-result-analysis.csv")
    ap.add_argument("oracles", help="directory of <model>-UBA.out vectors")
    ap.add_argument("-o", "--out")
    args = ap.parse_args()

    consensus = collections.defaultdict(dict)
    with open(args.raw, newline="") as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            if row[0] != "BVT-2026" or row[2] != "UpperBounds":
                continue
            for i, v in enumerate(ts.split_results(row[15])):
                if v not in ts.NO_ANSWER:
                    consensus[row[1]][i] = ts.normalize(v)

    census = collections.Counter()
    rows = []
    for name in sorted(os.listdir(args.oracles)):
        if not name.endswith("-UBA.out"):
            continue
        model = name[:-8]
        tgz = os.path.join(args.inputs, model + ".tgz")
        if not os.path.exists(tgz):
            census["no archive"] += 1
            continue
        places, single = read_archive(tgz, model)
        index = {p: i for i, p in enumerate(places)}
        vec = read_vector(os.path.join(args.oracles, name))
        if len(vec) != len(places):
            print(f"{model}: vector has {len(vec)} tokens, PNML {len(places)} places", file=sys.stderr)
            census["length mismatch"] += 1
            continue
        for fi, place in single.items():
            cons = consensus.get(model, {}).get(fi)
            if cons is None or place not in index:
                continue
            ours = vec[index[place]]
            if ours == "?":
                status = "vector undecided"
            elif cons == "+inf":
                status = "consensus inf, vector " + ours
            elif ours == cons:
                status = "confirmed"
            else:
                status = "CONTRADICTION"
                print(f"CONTRADICTION {model} formula {fi:02d} place {place} (p{index[place]}): vector {ours}, consensus {cons}")
            census[status] += 1
            rows.append({"Model": model, "formula": fi, "place": place, "index": index[place], "vector": ours, "consensus": cons, "status": status})
    for k, v in census.most_common():
        print(f"{k:30s} {v}", file=sys.stderr)
    if args.out:
        with open(args.out, "w", newline="") as f:
            w = csv.DictWriter(f, ["Model", "formula", "place", "index", "vector", "consensus", "status"])
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    main()
