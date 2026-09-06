#!/usr/bin/env python3
"""Turn the harness logs of the total examinations into one CSV row per log.

A total examination (QuasiLivenessAll, StableMarkingAll, UpperBoundsAll) asks
one question per transition or per place; the tool prints a `TOTAL <exam> <n>`
header then one line per atom as it is answered, `QLIVE t12 TRUE <techniques>`,
`BOUND p3 7 <techniques>`, and `BOUND p3 ? lo hi` while a bound is open. The
harness wraps each line in a teamcity marker whose duration is the time since
the previous one, which is how the closing curve is recovered.

    total-runs.csv   one row per log

With --oracles DIR, the answers are also written as vector oracles in the
pnmcc-models-2026 format (`<model>-QLA.out`...), a `?` where the run said
nothing, ready to be diffed against or merged into the published ones.
"""

import argparse
import csv
import os
import re
import sys

import mcclogs2csv as mcc

COLUMNS = [
    "log", "Model", "Examination", "host", "version", "timeout(s)", "start",
    "duration(ms)", "failure",
    "atoms", "answered", "completion", "witnessed", "proved", "open",
    "initial", "walk", "smt", "dd", "other",
    "t25(ms)", "t50(ms)", "t75(ms)", "t100(ms)",
    "petrispot", "walk calls", "walk asked", "walk solved", "walk ms",
    "last",
]

RE_HEADER = re.compile(r"^TOTAL (\w+) (\d+)")
RE_ATOM = re.compile(r"^(QLIVE|STABLE|BOUND) ([pt])(\d+) (\S+)(.*)$")

SUFFIX = {"QuasiLivenessAll": "QLA", "StableMarkingAll": "SMA", "UpperBoundsAll": "UBA"}
KEYWORD = {"QuasiLivenessAll": "QLIVE", "StableMarkingAll": "STABLE", "UpperBoundsAll": "BOUND"}

# The self certifying side of each examination: a firing sequence shows it.
WITNESSED = {"QLIVE": "TRUE", "STABLE": "FALSE"}


def technique(tags):
    """One family per verdict, the most specific engine that could have produced it."""
    if "DECISION_DIAGRAMS" in tags:
        return "dd"
    if "SAT_SMT" in tags:
        return "smt"
    if "WALK" in tags:
        return "walk"
    if "INITIAL_STATE" in tags:
        return "initial"
    return "other"


def parse(path):
    run = {c: "" for c in COLUMNS}
    run["log"] = path
    atoms = 0
    verdict = {}     # index -> (value, techniques)
    open_bounds = set()
    when = []        # cumulative ms at which each new verdict landed
    clock = 0
    petrispot = 0
    walk = [0, 0, 0, 0]
    failure = ""
    last = ""        # the tool's last line before the harness trailer
    trailer = False

    with open(path, errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            m = mcc.RE_SYSCALL.match(line)
            if m:
                words = m.group(1).split()
                if len(words) >= 3:
                    run["Model"], run["Examination"] = words[1], words[2]
                continue
            m = RE_HEADER.match(line)
            if m:
                atoms = int(m.group(2))
                continue
            m = RE_ATOM.match(line)
            if m:
                idx, value, tags = int(m.group(3)), m.group(4), m.group(5)
                if value == "?":
                    if idx not in verdict:
                        open_bounds.add(idx)
                elif idx not in verdict:
                    verdict[idx] = (value, tags)
                    open_bounds.discard(idx)
                    when.append(clock)
                continue
            if line.startswith("TIME LIMIT") or line.startswith("Actual read output values"):
                trailer = True
            if not trailer and line and not line.startswith("##teamcity") and not line.startswith(" Formula"):
                last = line[:70]
            m = mcc.RE_TC.search(line)
            if m:
                kind, name, duration = m.group(1), m.group(2), m.group(3)
                if kind == "testFinished" and duration is not None:
                    if name == "all":
                        run["duration(ms)"] = duration
                    else:
                        # runits is the time to the first line, the others the
                        # time since the previous line; the marker of a verdict
                        # line comes right after it, so the clock is read late
                        # by one line and corrected here
                        clock += int(duration)
                        if when and when[-1] == clock - int(duration) and name != "runits":
                            when[-1] = clock
                continue
            if line.startswith("Running PetriSpot"):
                petrispot += 1
                continue
            m = mcc.RE_WALK.match(line)
            if m:
                walk[0] += 1
                walk[1] += int(m.group(2))
                walk[2] += int(m.group(1))
                walk[3] += int(m.group(3))
                continue
            if not failure:
                for name, pattern in mcc.FAILURES:
                    if pattern in line:
                        failure = name
                        break
            if not run["start"]:
                m = mcc.RE_STAMP.match(line)
                if m:
                    run["start"] = m.group(1)
            for key, rx in (("host", mcc.RE_HOST), ("timeout(s)", mcc.RE_TIMEOUT),
                            ("version", mcc.RE_VERSION)):
                m = rx.match(line)
                if m:
                    run[key] = m.group(1)
                    break

    kw = KEYWORD.get(run["Examination"], "")
    census = {"initial": 0, "walk": 0, "smt": 0, "dd": 0, "other": 0}
    witnessed = 0
    for value, tags in verdict.values():
        census[technique(tags)] += 1
        if WITNESSED.get(kw) == value:
            witnessed += 1
    answered = len(verdict)
    run.update(census)
    run.update({
        "failure": failure, "atoms": atoms, "answered": answered,
        "completion": f"{answered / atoms:.4f}" if atoms else "",
        "witnessed": witnessed, "proved": answered - witnessed if kw != "BOUND" else "",
        "open": len(open_bounds),
        "petrispot": petrispot, "walk calls": walk[0], "walk asked": walk[1],
        "walk solved": walk[2], "walk ms": walk[3], "last": last,
    })
    for q, col in ((0.25, "t25(ms)"), (0.5, "t50(ms)"), (0.75, "t75(ms)"), (1.0, "t100(ms)")):
        need = int(atoms * q + 0.999999) if atoms else 0
        if need and answered >= need:
            run[col] = when[need - 1]
    return run, atoms, verdict


def write_oracle_to(f, run, atoms, verdict):
    """The answers as a vector oracle, wrapped at 80 columns, on an open file."""
    exam = run["Examination"]
    f.write(f"{run['Model']} {exam}\n{KEYWORD[exam]}\n")
    if exam == "UpperBoundsAll":
        toks = [verdict[i][0] if i in verdict else "?" for i in range(atoms)]
        line = ""
        for t in toks:
            if len(line) + len(t) + 1 > 80:
                f.write(line + "\n")
                line = t
            else:
                line = t if not line else line + " " + t
        f.write(line + "\n")
    else:
        chars = "".join(verdict[i][0][0] if i in verdict else "?" for i in range(atoms))
        for i in range(0, atoms, 80):
            f.write(chars[i:i + 80] + "\n")


def write_oracle(outdir, run, atoms, verdict):
    """The vector oracle as <model>-<QLA|SMA|UBA>.out in outdir."""
    exam = run["Examination"]
    if exam not in SUFFIX or not atoms:
        return
    path = os.path.join(outdir, f"{run['Model']}-{SUFFIX[exam]}.out")
    with open(path, "w") as f:
        write_oracle_to(f, run, atoms, verdict)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("dirs", nargs="+", help="directories holding *out logs, or logs")
    ap.add_argument("-o", "--outdir", default=".", help="where to write the CSV")
    ap.add_argument("--oracles", help="also write one vector oracle per log here")
    args = ap.parse_args()

    logs = []
    for d in args.dirs:
        if os.path.isdir(d):
            logs += sorted(os.path.join(d, f) for f in os.listdir(d) if f.endswith("out") or f.endswith("log"))
        else:
            logs.append(d)

    os.makedirs(args.outdir, exist_ok=True)
    if args.oracles:
        os.makedirs(args.oracles, exist_ok=True)
    runs_path = os.path.join(args.outdir, "total-runs.csv")
    with open(runs_path, "w", newline="") as rf:
        w = csv.DictWriter(rf, COLUMNS)
        w.writeheader()
        for path in logs:
            run, atoms, verdict = parse(path)
            w.writerow(run)
            if args.oracles:
                write_oracle(args.oracles, run, atoms, verdict)
    print(f"{len(logs)} logs -> {runs_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
