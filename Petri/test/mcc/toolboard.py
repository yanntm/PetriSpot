#!/usr/bin/env python3
"""Our campaign against every tool of the contest, per examination.

`raw-result-analysis.csv` of `pnmcc-models-20xx` holds one row per (tool,
model instance, examination) with the tool's verdicts; the oracle in our logs
is the consensus of that same file. Scoring each tool's tokens against the
consensus the way `mcclogs2csv.py` scores ours gives one table per
examination: values answered, of which correct, wrong, and bonus (the
consensus has no value). Our row comes from `verdicts.csv`.

With --tool, the instances where that tool got more values right than we did
are listed, with both counts and both wall times (ours from `runs.csv` next to
the verdicts), which is where to look for what the contest build did better.
Contest runs had 3600 s and 4 cores; ours are whatever `runs.csv` says.

    toolboard.py raw-result-analysis.csv verdicts.csv [--tool ITS-Tools] [-o toolboard.csv]
"""

import argparse
import collections
import csv
import os
import sys

import toolsupport as ts


def read_ours(path):
    """(model, exam) -> {index: (oracle, answer, status)} from verdicts.csv."""
    ours = collections.defaultdict(dict)
    exams = []
    with open(path, newline="") as f:
        for v in csv.DictReader(f):
            key = (v["Model"], v["Examination"])
            if v["Examination"] not in exams:
                exams.append(v["Examination"])
            ours[key][ts.index_of(v["formula"], v["Examination"])] = (v["oracle"], v["answer"], v["status"])
    return ours, exams


def read_tools(path, exams):
    """(tool, model, exam) -> (list of normalized tokens, time ms)."""
    tools = {}
    with open(path, newline="") as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            tool, model, exam, results, time = row[0], row[1], row[2], row[6], row[10]
            if exam not in exams:
                continue
            toks = [ts.normalize(v) if v not in ts.NO_ANSWER else "" for v in ts.split_results(results)]
            tools[(tool, model, exam)] = (toks, time)
    return tools


def score(tokens, oracle):
    """answered, correct, wrong, bonus of a token list against an oracle {index: value}."""
    a = c = w = b = 0
    for i, t in enumerate(tokens):
        if not t:
            continue
        a += 1
        o = oracle.get(i, "?")
        if o == "?":
            b += 1
        elif t == o:
            c += 1
        else:
            w += 1
    return a, c, w, b


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("raw", help="raw-result-analysis.csv")
    ap.add_argument("verdicts", help="verdicts.csv from mcclogs2csv.py")
    ap.add_argument("--tool", help="list the instances where this tool beat us")
    ap.add_argument("-o", "--out", help="write the per tool table here as CSV")
    args = ap.parse_args()

    ours, exams = read_ours(args.verdicts)
    tools = read_tools(args.raw, exams)
    runs = {}
    runs_path = os.path.join(os.path.dirname(args.verdicts), "runs.csv")
    if os.path.exists(runs_path):
        with open(runs_path, newline="") as f:
            for r in csv.DictReader(f):
                runs[(r["Model"], r["Examination"])] = r

    names = sorted({t for t, _, _ in tools} - ts.VIRTUAL)
    table = []
    for exam in exams:
        instances = sorted(m for m, e in ours if e == exam)
        rows = []
        for tool in names:
            tot = [0, 0, 0, 0]
            for m in instances:
                oracle = {i: o for i, (o, _, _) in ours[(m, exam)].items()}
                toks, _ = tools.get((tool, m, exam), ([], ""))
                for k, v in enumerate(score(toks, oracle)):
                    tot[k] += v
            rows.append((tool, *tot))
        mine = [0, 0, 0, 0]
        for m in instances:
            for _, (o, a, status) in ours[(m, exam)].items():
                if status in ("ok", "wrong", "bonus"):
                    mine[0] += 1
                mine[1] += status == "ok"
                mine[2] += status == "wrong"
                mine[3] += status == "bonus"
        rows.append(("us", *mine))
        print(f"\n## {exam} ({len(instances)} instances)")
        print(f"{'tool':12s} {'answered':>9s} {'correct':>8s} {'wrong':>6s} {'bonus':>6s}")
        for tool, a, c, w, b in sorted(rows, key=lambda r: -r[2]):
            print(f"{tool:12s} {a:9d} {c:8d} {w:6d} {b:6d}")
            table.append({"Examination": exam, "tool": tool, "answered": a, "correct": c, "wrong": w, "bonus": b})

        if args.tool:
            print(f"\n   instances where {args.tool} got more right than we did (their/our correct, their/our time s):")
            n = 0
            for m in instances:
                oracle = {i: o for i, (o, _, _) in ours[(m, exam)].items()}
                toks, time = tools.get((args.tool, m, exam), ([], ""))
                theirs = score(toks, oracle)[1]
                us = sum(1 for _, (_, _, s) in ours[(m, exam)].items() if s == "ok")
                if theirs > us:
                    n += 1
                    r = runs.get((m, exam), {})
                    dur = int(r.get("duration(ms)") or 0) / 1000
                    fail = r.get("failure", "")
                    print(f"   {m:42s} {theirs:3d}/{us:<3d} {float(time or 0)/1000:8.1f}/{dur:<8.1f} {fail}")
            print(f"   {n} instances")

    if args.out:
        with open(args.out, "w", newline="") as f:
            w = csv.DictWriter(f, ["Examination", "tool", "answered", "correct", "wrong", "bonus"])
            w.writeheader()
            w.writerows(table)
        print(f"{len(table)} rows -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
