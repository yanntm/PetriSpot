"""The report sections of the total examinations, from `total-runs.csv`.

One section per examination present: headline, completion deciles, time
shape, the families with incomplete runs, the lowest completions at the wall,
the wall runs that never reached a walk, failures, and the completion changes
against a baseline `total-runs.csv` when one is given. Used by `report.py`.
"""

import collections
import csv
import statistics

from mdtable import I, family, table, hours, seconds, at_wall, bullet

EXAMS = ("QuasiLivenessAll", "StableMarkingAll", "UpperBoundsAll")
ENGINES = ("initial", "walk", "smt", "dd", "other")


def read(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def completion(r):
    return float(r["completion"]) if r["completion"] else 0.0


def section(rows, exam, timeout_s, baseline=None):
    rs = [r for r in rows if r["Examination"] == exam]
    if not rs:
        return ""
    live = [r for r in rs if I(r["atoms"])]
    out = [f"## {exam}\n"]
    atoms = sum(I(r["atoms"]) for r in rs)
    answered = sum(I(r["answered"]) for r in rs)
    complete = [r for r in live if completion(r) >= 1]
    incomplete = [r for r in live if completion(r) < 1]
    wall = [r for r in rs if at_wall(r, timeout_s)]
    wit = sum(I(r["witnessed"]) for r in rs)
    prov = sum(I(r["proved"]) for r in rs)
    opn = sum(I(r["open"]) for r in rs)
    eng = collections.Counter()
    for r in rs:
        for k in ENGINES:
            eng[k] += I(r[k])
    hl = [("runs", len(rs)), ("atoms", atoms), ("answered", answered),
          (f"completion", f"{answered / atoms:.3f}" if atoms else ""),
          ("runs complete", len(complete)), ("runs incomplete", len(incomplete)),
          ("runs at the wall", len(wall)),
          ("wall hours", f"{sum(hours(r['duration(ms)']) for r in rs):.0f}"),
          ("walker hours", f"{sum(hours(r['walk ms']) for r in rs):.1f}"),
          ("walker calls", sum(I(r["walk calls"]) for r in rs))]
    if exam == "UpperBoundsAll":
        hl.append(("bounds still open", opn))
    else:
        hl += [("witnessed", wit), ("proved", prov)]
    out.append(table(["", ""], hl, "lr"))
    out.append("\nVerdicts by engine: " + ", ".join(f"{k} {eng[k]}" for k in ENGINES) + ".\n")

    # deciles
    hist = collections.Counter(min(int(completion(r) * 10), 9) for r in live)
    out.append("\n### Completion of the runs\n\n")
    out.append(table(["completion", "runs"], [(f"{d / 10:.1f} .. {(d + 1) / 10:.1f}" if d < 9 else "0.9 .. 1.0", hist.get(d, 0)) for d in range(10)], "lr"))

    # time shape
    out.append("\n### Time shape\n\n")
    if complete:
        d = sorted(seconds(r["duration(ms)"]) for r in complete)
        t100 = [seconds(r["t100(ms)"]) for r in complete if r["t100(ms)"]]
        out.append(f"Complete runs: duration median {statistics.median(d):.0f} s, p90 {d[int(.9 * (len(d) - 1))]:.0f} s, max {max(d):.0f} s; "
                   f"last atom closed at median {statistics.median(t100):.1f} s.\n\n")
    wq = collections.Counter(sum(1 for c in ("t25(ms)", "t50(ms)", "t75(ms)") if r[c]) for r in wall)
    out.append("Runs at the wall by quartiles reached: " + ", ".join(f"{wq.get(k, 0)} reached {k}/3" for k in range(4)) + ".\n\n")
    t75 = sorted(seconds(r["t75(ms)"]) for r in wall if r["t75(ms)"])
    if t75:
        out.append(f"Among wall runs that closed three quarters: that point came at median {statistics.median(t75):.0f} s "
                   f"(p25 {t75[len(t75) // 4]:.0f} s, p75 {t75[3 * len(t75) // 4]:.0f} s): a knee, then a residue the rest of the time does not close.\n\n")
    early = [r for r in incomplete if not at_wall(r, timeout_s)]
    if early:
        fails = collections.Counter(r["failure"] or "ended early" for r in early)
        out.append(f"Incomplete runs that did not reach the wall: {len(early)} ({', '.join(f'{v} {k}' for k, v in fails.most_common())}).\n\n")

    # families
    byf = collections.defaultdict(list)
    for r in live:
        byf[family(r["Model"])].append(r)
    fam_rows = []
    for f, L in byf.items():
        comp = [completion(r) for r in L]
        if min(comp) < 1:
            una = sum(I(r["atoms"]) - I(r["answered"]) for r in L)
            fam_rows.append((f, len(L), sum(1 for c in comp if c >= 1), sum(comp) / len(comp), una, sum(hours(r["duration(ms)"]) for r in L)))
    fam_rows.sort(key=lambda x: -x[4])
    out.append(f"### Families with incomplete runs ({len(fam_rows)} of {len(byf)})\n\n")
    out.append(table(["family", "instances", "complete", "mean completion", "unanswered atoms", "wall h"], fam_rows, "lrrrrr"))

    # lowest completion at the wall
    out.append("\n### Lowest completion among runs at the wall\n\n")
    low = sorted(wall, key=completion)[:30]
    out.append(table(["instance", "completion", "atoms", "answered", "initial", "walk", "smt", "dd", "walk calls", "walk s", "last activity"],
                     [(r["Model"], completion(r), I(r["atoms"]), I(r["answered"]), I(r["initial"]), I(r["walk"]), I(r["smt"]), I(r["dd"]),
                       I(r["walk calls"]), int(seconds(r["walk ms"])), r.get("last", "")[:60].replace("|", "/")) for r in low],
                     "lrrrrrrrrrl"))

    # wall runs with no walk
    nowalk = [r for r in wall if I(r["walk calls"]) == 0]
    out.append(f"\n### Runs at the wall that never called the walker: {len(nowalk)}\n\n")
    if nowalk:
        lastc = collections.Counter(r.get("last", "")[:45] for r in nowalk)
        out.append("By last activity line:\n\n")
        out.append(table(["last activity", "runs"], [(k.replace("|", "/"), v) for k, v in lastc.most_common(10)], "lr"))
        out.append("\n" + bullet([f"{r['Model']} ({I(r['answered'])}/{I(r['atoms'])})" for r in sorted(nowalk, key=lambda r: r["Model"])], 60))

    # failures
    fails = collections.Counter(r["failure"] for r in rs if r["failure"])
    if fails:
        out.append("\n### Failures\n\n")
        out.append(table(["failure", "runs"], fails.most_common(), "lr"))
        out.append("\n" + bullet([f"{r['Model']}: {r['failure']}" for r in rs if r["failure"]]))

    # baseline
    if baseline:
        base = {r["Model"]: r for r in baseline if r["Examination"] == exam}
        changes = []
        for r in live:
            b = base.get(r["Model"])
            if b and abs(completion(b) - completion(r)) >= 0.01:
                changes.append((r["Model"], completion(b), completion(r), I(b["answered"]), I(r["answered"]), int(seconds(b["duration(ms)"])), int(seconds(r["duration(ms)"]))))
        changes.sort(key=lambda x: x[1] - x[2])
        ans_b = sum(I(b["answered"]) for b in base.values())
        out.append(f"\n### Against the baseline\n\nAnswered {answered} atoms against {ans_b} ({answered - ans_b:+d}). "
                   f"{len(changes)} instances moved by at least a percent:\n\n")
        out.append(table(["instance", "baseline", "now", "answered before", "after", "s before", "s after"], changes[:60], "lrrrrrr"))
    return "".join(out) + "\n"


def sections(path, timeout_s, baseline_path=None):
    rows = read(path)
    base = read(baseline_path) if baseline_path else None
    return "".join(section(rows, ex, timeout_s, base) for ex in EXAMS)
