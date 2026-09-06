#!/usr/bin/env python3
"""A markdown report of a collected campaign, to read and to argue from.

Takes a folder written by the collectors (`runs.csv`, `verdicts.csv`,
`support.csv` from `mcclogs2csv.py` and `toolsupport.py`; `total-runs.csv`
from `totallogs2csv.py`), and writes `REPORT.md` there: one section per
examination with the headline census, the wrong verdicts, the missed values
and who else has them, the bonus values, the runs at the wall, the runs that
gave up early, the failures; against the contest tools when the raw results
are given; against another campaign folder when a baseline is given. The
total examinations get their own sections (`report_total.py`).

    report.py csv/2026-09-06/ [--baseline csv/2026-09-05/] [--raw raw-result-analysis.csv] [--title ...]

Numbers come from the CSVs only; rerun the collectors first when the logs
changed. `pandoc REPORT.md -o REPORT.pdf` if paper is wanted.
"""

import argparse
import collections
import csv
import os
import sys

import report_total
from mdtable import I, family, table, hours, seconds, at_wall, truncated, bullet

STATUS = ("ok", "wrong", "missed", "none", "bonus", "extra")


def read(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def timeout_of(runs):
    t = collections.Counter(r["timeout(s)"] for r in runs if r["timeout(s)"])
    return int(t.most_common(1)[0][0]) if t else 1800


def headline(runs, timeout_s):
    c = collections.Counter()
    for r in runs:
        for k in ("known", "answered") + STATUS:
            c[k] += I(r[k])
    rows = [("instances", len(runs)),
            ("build", ", ".join(f"{v} ({n})" for v, n in collections.Counter(r["version"] for r in runs).most_common())),
            ("timeout s", timeout_s),
            ("wall hours", f"{sum(hours(r['duration(ms)']) for r in runs):.1f}"),
            ("runs at the wall", sum(1 for r in runs if at_wall(r, timeout_s))),
            ("truncated logs (no trailer)", sum(1 for r in runs if truncated(r))),
            ("runs with a failure signature", sum(1 for r in runs if r["failure"])),
            ("oracle values known", c["known"]), ("answered", c["answered"]),
            ("ok", c["ok"]), ("**wrong**", f"**{c['wrong']}**"), ("missed (the field has it)", c["missed"]),
            ("none (nobody has it)", c["none"]), ("bonus (only we have it)", c["bonus"]),
            ("walker calls / solved", f"{sum(I(r['walk calls']) for r in runs)} / {sum(I(r['walk solved']) for r in runs)}"),
            ("walker hours (as reported by the calls)", f"{sum(hours(r['walk ms']) for r in runs):.1f}")]
    return table(["", ""], rows, "lr")


def section(exam, runs, verdicts, support, timeout_s, baseline, raw):
    out = [f"## {exam}\n\n", headline(runs, timeout_s), "\n"]
    sup = {(s["Model"], s["formula"]): s for s in support}
    byrun = {r["Model"]: r for r in runs}

    wrong = [v for v in verdicts if v["status"] == "wrong"]
    out.append(f"### Wrong verdicts: {len(wrong)}\n\n")
    if wrong:
        out.append(table(["formula", "oracle", "ours", "backed by", "run s"],
                         [(v["formula"], v["oracle"], v["answer"], sup.get((v["Model"], v["formula"]), {}).get("who", ""),
                           int(seconds(byrun[v["Model"]]["duration(ms)"]))) for v in wrong], "lllll"))

    missed = [v for v in verdicts if v["status"] == "missed"]
    out.append(f"\n### Missed values: {len(missed)}\n\n")
    if missed:
        who = collections.Counter()
        single = collections.Counter()
        for v in missed:
            s = sup.get((v["Model"], v["formula"]), {})
            ts = s.get("who", "").split()
            for t in ts:
                who[t] += 1
            if len(ts) == 1:
                single[ts[0]] += 1
        out.append("Backed by: " + ", ".join(f"{t} {n}" for t, n in who.most_common()) + ". ")
        out.append("Backed by a single tool: " + (", ".join(f"{t} {n}" for t, n in single.most_common()) or "none") + ".\n\n")
        kinds = collections.Counter()
        for v in missed:
            r = byrun[v["Model"]]
            kinds["in a run with a failure signature" if r["failure"] else "in a truncated log" if truncated(r)
                  else "in a run at the wall" if at_wall(r, timeout_s) else "in a run that ended early"] += 1
        out.append("Where: " + ", ".join(f"{n} {k}" for k, n in kinds.most_common()) + ".\n\n")
        fam = collections.Counter(family(v["Model"]) for v in missed)
        out.append(table(["family", "missed"], fam.most_common(20), "lr"))
        out.append("\nMissed values that only the contest ITS-Tools produced (our own regressions):\n\n")
        out.append(bullet([f"{v['formula']} = {v['oracle']}, run {int(seconds(byrun[v['Model']]['duration(ms)']))} s {byrun[v['Model']]['failure']}"
                           for v in missed if sup.get((v["Model"], v["formula"]), {}).get("who", "").strip() == "ITS-Tools"], 40) or "none\n")

    bonus = [v for v in verdicts if v["status"] == "bonus"]
    out.append(f"\n### Bonus values (the consensus has none): {len(bonus)}\n\n")
    out.append(bullet([f"{v['formula']} = {v['answer']}" for v in bonus], 40))

    wall = [r for r in runs if at_wall(r, timeout_s)]
    out.append(f"\n### Runs at the wall: {len(wall)}\n\n")
    out.append(table(["instance", "ok", "missed", "none", "walk calls", "walk s"],
                     [(r["Model"], I(r["ok"]), I(r["missed"]), I(r["none"]), I(r["walk calls"]), int(seconds(r["walk ms"])))
                      for r in sorted(wall, key=lambda r: (-I(r["missed"]), r["Model"]))[:80]], "lrrrrr"))
    if len(wall) > 80:
        out.append(f"\n... {len(wall) - 80} more\n")

    early = [r for r in runs if (I(r["missed"]) + I(r["none"])) > 0 and not at_wall(r, timeout_s) and not r["failure"] and not truncated(r)]
    out.append(f"\n### Runs that ended early with open values and no failure signature: {len(early)}\n\n")
    out.append(table(["instance", "s", "ok", "missed", "none", "walk calls / solved"],
                     [(r["Model"], int(seconds(r["duration(ms)"])), I(r["ok"]), I(r["missed"]), I(r["none"]), f"{r['walk calls']}/{r['walk solved']}")
                      for r in sorted(early, key=lambda r: r["Model"])[:60]], "lrrrrl"))

    fails = [r for r in runs if r["failure"] or truncated(r)]
    out.append(f"\n### Failures and truncated logs: {len(fails)}\n\n")
    out.append(table(["instance", "failure", "s", "ok", "missed"],
                     [(r["Model"], r["failure"] or "truncated log", int(seconds(r["duration(ms)"])), I(r["ok"]), I(r["missed"])) for r in sorted(fails, key=lambda r: r["Model"])], "llrrr"))

    if raw:
        out.append(board(exam, verdicts, runs, raw))
    if baseline:
        out.append(against(exam, verdicts, runs, baseline))
    return "".join(out) + "\n"


def board(exam, verdicts, runs, raw):
    """Every contest tool scored against the consensus on this examination, and where ITS-Tools 2026 beat us."""
    import toolboard as tb
    import toolsupport as ts
    ours = collections.defaultdict(dict)
    for v in verdicts:
        ours[v["Model"]][ts.index_of(v["formula"], exam)] = (v["oracle"], v["answer"], v["status"])
    tools = tb.read_tools(raw, [exam])
    names = sorted({t for t, _, e in tools if e == exam} - ts.VIRTUAL)
    rows = []
    for tool in names:
        tot = [0, 0, 0, 0]
        for m in ours:
            oracle = {i: o for i, (o, _, _) in ours[m].items()}
            for k, v in enumerate(tb.score(tools.get((tool, m, exam), ([], ""))[0], oracle)):
                tot[k] += v
        rows.append((tool, *tot))
    mine = [sum(1 for m in ours for _, (_, _, s) in ours[m].items() if s in ("ok", "wrong", "bonus"))] + \
           [sum(1 for m in ours for _, (_, _, s) in ours[m].items() if s == st) for st in ("ok", "wrong", "bonus")]
    rows.append(("**us**", *mine))
    rows.sort(key=lambda r: -r[2])
    out = [f"\n### Against the contest tools (they had 3600 s)\n\n", table(["tool", "answered", "correct", "wrong", "bonus"], rows, "lrrrr")]
    byrun = {r["Model"]: r for r in runs}
    beat = []
    for m in sorted(ours):
        oracle = {i: o for i, (o, _, _) in ours[m].items()}
        toks, time = tools.get(("ITS-Tools", m, exam), ([], ""))
        theirs = tb.score(toks, oracle)[1]
        us = sum(1 for _, (_, _, s) in ours[m].items() if s == "ok")
        if theirs > us:
            r = byrun.get(m, {})
            beat.append((m, f"{theirs}/{us}", int(float(time or 0) / 1000), int(seconds(r.get("duration(ms)"))), r.get("failure", "")))
    out.append(f"\nInstances where the contest ITS-Tools got more right than we did: {len(beat)}\n\n")
    out.append(table(["instance", "their/our correct", "their s", "our s", "our failure"], beat, "llrrl"))
    return "".join(out)


def against(exam, verdicts, runs, baseline):
    """Per formula status changes against another campaign folder."""
    bv = {(v["Model"], v["formula"]): v for v in read(os.path.join(baseline, "verdicts.csv")) if v["Examination"] == exam}
    br = {r["Model"]: r for r in read(os.path.join(baseline, "runs.csv")) if r["Examination"] == exam}
    if not bv:
        return ""
    byrun = {r["Model"]: r for r in runs}
    changes = []
    for v in verdicts:
        b = bv.get((v["Model"], v["formula"]))
        if b and (b["status"] != v["status"] or b["answer"] != v["answer"]):
            changes.append((v["Model"], v["formula"].replace(v["Model"] + "-", ""), b["status"], b["answer"], int(seconds(br[v["Model"]]["duration(ms)"])) if v["Model"] in br else "",
                            v["status"], v["answer"], int(seconds(byrun[v["Model"]]["duration(ms)"])), byrun[v["Model"]]["failure"]))
    rank = {"wrong": 0, "missed": 1, "none": 2, "bonus": 3, "ok": 4}
    changes.sort(key=lambda c: (rank.get(c[5], 5) - rank.get(c[2], 5), c[0]))
    census = collections.Counter(f"{c[2]} -> {c[5]}" for c in changes)
    out = [f"\n### Against the baseline `{baseline}`\n\n", "Status changes: " + (", ".join(f"{k} {n}" for k, n in census.most_common()) or "none") + ".\n\n"]
    out.append(table(["instance", "formula", "was", "answer", "s", "now", "answer", "s", "failure"], changes[:100], "lllllllll"))
    if len(changes) > 100:
        out.append(f"\n... {len(changes) - 100} more\n")
    return "".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("folder", help="a folder written by the collectors")
    ap.add_argument("--baseline", help="another such folder to compare against")
    ap.add_argument("--raw", help="raw-result-analysis.csv of the contest, for the tool board")
    ap.add_argument("--title", default=None)
    ap.add_argument("-o", "--out", help="default <folder>/REPORT.md")
    args = ap.parse_args()

    runs = read(os.path.join(args.folder, "runs.csv"))
    verdicts = read(os.path.join(args.folder, "verdicts.csv"))
    support = read(os.path.join(args.folder, "support.csv"))
    timeout_s = timeout_of(runs) if runs else 1800
    title = args.title or f"Campaign report: {os.path.basename(os.path.normpath(args.folder))}"
    out = [f"# {title}\n\n", f"Generated by `report.py` from `{args.folder}`"]
    if args.baseline:
        out.append(f", against `{args.baseline}`")
    out.append(". Vocabulary: `missed` is a value the consensus has and we do not, `none` a value nobody has, `bonus` a value only we have; "
               "a run is *at the wall* when the harness killed it at the timeout.\n\n")
    exams = []
    for r in runs:
        if r["Examination"] not in exams:
            exams.append(r["Examination"])
    for exam in exams:
        out.append(section(exam, [r for r in runs if r["Examination"] == exam], [v for v in verdicts if v["Examination"] == exam],
                           [s for s in support if s["Examination"] == exam], timeout_s, args.baseline, args.raw))
    total = os.path.join(args.folder, "total-runs.csv")
    if os.path.exists(total):
        trs = report_total.read(total)
        t = timeout_of(trs) if trs else timeout_s
        bt = os.path.join(args.baseline, "total-runs.csv") if args.baseline else None
        out.append(report_total.sections(total, t, bt if bt and os.path.exists(bt) else None))
    path = args.out or os.path.join(args.folder, "REPORT.md")
    with open(path, "w") as f:
        f.write("".join(out))
    print(f"{path}: {sum(len(s) for s in out)} chars, examinations {exams}", file=sys.stderr)


if __name__ == "__main__":
    main()
