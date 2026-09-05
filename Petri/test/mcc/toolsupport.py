#!/usr/bin/env python3
"""Who else answered the values we did not.

`raw-result-analysis.csv` of `pnmcc-models-20xx` holds one row per (tool, model
instance, examination), with the tool's verdicts, the consensus verdicts and a
count of the tools that back the consensus. The verdicts are positional: the
i-th token of the row answers the i-th formula of the examination, which is how
`csv_to_control.pl` fills the oracle skeletons.

Joining that on a `verdicts.csv` produced by `mcclogs2csv.py` says, for every
value we did not produce, which tools did.

    toolsupport.py raw-result-analysis.csv verdicts.csv [-o support.csv]
"""

import argparse
import collections
import csv
import re
import sys

# `BVT-2026`, the best virtual tool, is the union of the others and so never
# adds evidence of its own. `2025-gold` is the previous edition's gold medal
# rerun here: a real tool, and for a few values the only one that answered.
VIRTUAL = {"BVT-2026"}

# The spellings of "no verdict": did not compute, did not finish, cannot
# compute, and the unanswered marker of the consensus column.
NO_ANSWER = {"DNC", "DNF", "CC", "?", ""}

RE_INDEX = re.compile(r"-(\d\d)$")


def normalize(v):
    """The verdict spellings of the raw results, as the oracle writes them."""
    v = v.strip()
    if v == "T":
        return "TRUE"
    if v == "F":
        return "FALSE"
    v = re.sub(r"(\d)\.0000E\+0005", r"\g<1>00000", v)
    if v.endswith("inf"):
        return "+inf"
    return v


def split_results(field):
    """The per-formula tokens of a `results` or `estimated result` field."""
    field = field.strip().replace("(", "").replace(")", "")
    if not field:
        return []
    if " " in field:
        return field.split()
    # A 16 character string is one character per formula, anything else is one
    # verdict for the whole examination (the global properties).
    if len(field) == 16:
        return list(field)
    return [field]


def read_raw(path):
    """(model, examination, index) -> {verdict: [tools that answered it]}"""
    support = collections.defaultdict(lambda: collections.defaultdict(list))
    with open(path, newline="") as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            tool, model, exam, results = row[0], row[1], row[2], row[6]
            if tool in VIRTUAL:
                continue
            for i, v in enumerate(split_results(results)):
                if v in NO_ANSWER:
                    continue
                support[(model, exam, i)][normalize(v)].append(tool)
    return support


def index_of(formula, examination):
    """The position of a formula in its examination, as the raw results order it."""
    m = RE_INDEX.search(formula)
    return int(m.group(1)) if m else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("raw", help="raw-result-analysis.csv")
    ap.add_argument("verdicts", help="verdicts.csv from mcclogs2csv.py")
    ap.add_argument("-o", "--out", default="support.csv")
    args = ap.parse_args()

    support = read_raw(args.raw)
    counts = collections.Counter()
    rows = []
    with open(args.verdicts, newline="") as f:
        for v in csv.DictReader(f):
            key = (v["Model"], v["Examination"], index_of(v["formula"], v["Examination"]))
            by_verdict = support.get(key, {})
            tools = sorted({t for ts in by_verdict.values() for t in ts})
            oracle = v["oracle"]
            backing = sorted(by_verdict.get(oracle, []))
            dissent = sorted(t for value, ts in by_verdict.items()
                             if value != oracle for t in ts)
            counts[(v["status"], len(backing))] += 1
            rows.append({
                "Model": v["Model"], "Examination": v["Examination"],
                "formula": v["formula"], "oracle": oracle,
                "answer": v["answer"], "status": v["status"],
                "tools": len(tools), "backing": len(backing),
                "who": " ".join(backing), "dissent": " ".join(dissent),
            })

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, ["Model", "Examination", "formula", "oracle",
                               "answer", "status", "tools", "backing", "who",
                               "dissent"])
        w.writeheader()
        w.writerows(rows)

    print(f"{len(rows)} values -> {args.out}", file=sys.stderr)
    for (status, n), c in sorted(counts.items()):
        print(f"{status:8s} backed by {n} tools : {c}", file=sys.stderr)


if __name__ == "__main__":
    main()
