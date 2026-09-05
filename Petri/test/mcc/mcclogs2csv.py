#!/usr/bin/env python3
"""Turn a directory of MCC harness logs into two CSV tables.

One log is one `run_test.pl` invocation: one model instance, one examination.
It carries the oracle verdicts (the `Control values :` block), the tool output
(`FORMULA` / `STATE_SPACE` lines) and the teamcity markers the harness prints.

    runs.csv     one row per log, with the counts of the verdict comparison
    verdicts.csv one row per (model, examination, formula)

The comparison mirrors `run_test.pl`: the last verdict a formula receives is
the one that counts, an oracle `?` is not a constraint, and `inf` in the oracle
matches `+inf` in the output.
"""

import argparse
import csv
import os
import re
import sys

RUN_COLUMNS = [
    "log", "Model", "Examination", "host", "version", "timeout(s)", "args",
    "start", "Test started", "Test fail", "Test fin", "duration(ms)",
    "first(ms)", "petrispot", "petrispot fail", "petrispot timeout",
    "walk calls", "walk asked", "walk solved", "walk ms", "walk killed",
    "failure",
    "oracle", "known", "answered", "ok", "wrong", "missed", "bonus", "extra",
]

VERDICT_COLUMNS = [
    "log", "Model", "Examination", "formula", "oracle", "answer", "status",
    "techniques",
]

RE_HOST = re.compile(r"^Linux (\S+) ")
RE_TITLE = re.compile(r"^Running test : (\S+)")
RE_TIMEOUT = re.compile(r"^Timeout set at :(\d+) seconds")
RE_SYSCALL = re.compile(r"^syscalling : (.*)$")
RE_VERSION = re.compile(r"^Running Version (\S+)")
RE_VERDICT = re.compile(r"^(?:FORMULA|STATE_SPACE)\s+(\S+)\s+(\S+)(?:\s+TECHNIQUES\s+(.*?))?\s*$")
RE_TC = re.compile(r"##teamcity\[(\w+) name='([^']*)'(?:.*?duration='(\d+)')?")
RE_STAMP = re.compile(r"^\[(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)\]")
RE_WALK = re.compile(r"^PetriSpot walker: (\d+)/(\d+) properties solved"
                     r"(?:, \d+ bounds reported)? in (\d+) ms \(exit (-?\d+)\)")

# What made a run produce nothing. The first pattern that matches names the
# failure; they are ordered from the most specific to the most generic.
FAILURES = [
    ("no_input", "Cannot open file"),
    ("overlarge_marking", "OverlargeMarkingException"),
    ("eclipse_fatal", "An error has occurred. See the log file"),
    ("out_of_memory", "OutOfMemoryError"),
    ("its_abort", "terminate called"),
]


def equivalent(expected, got):
    """The oracle accepts `got` for `expected`, as `run_test.pl` decides it."""
    if expected == "?":
        return True
    if expected == "inf" and got == "+inf":
        return True
    return expected == got


def parse(path):
    """Read one log; return its run row and its list of verdict rows."""
    run = {c: "" for c in RUN_COLUMNS}
    run["log"] = path
    started = failed = finished = 0
    oracle = {}
    answer = {}
    techniques = {}
    in_control = False
    seen_syscall = False
    petrispot = petrispot_fail = petrispot_timeout = 0
    walk = [0, 0, 0, 0, 0]  # calls, asked, solved, ms, killed
    failure = ""

    with open(path, errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if in_control:
                if "=" in line and not line.startswith("syscalling"):
                    key, _, value = line.partition("=")
                    oracle[key] = value
                    continue
                in_control = False
            if line.startswith("Control values :"):
                in_control = True
                continue
            m = RE_SYSCALL.match(line)
            if m:
                seen_syscall = True
                words = m.group(1).split()
                # ./runatest[_cluster].sh <model> <examination> [args...]
                if len(words) >= 3:
                    run["Model"] = words[1]
                    run["Examination"] = words[2]
                    run["args"] = " ".join(words[3:])
                continue
            if seen_syscall:
                m = RE_VERDICT.match(line)
                if m:
                    answer[m.group(1)] = m.group(2)
                    techniques[m.group(1)] = (m.group(3) or "").strip()
                    continue
            if line.startswith("Running PetriSpot"):
                petrispot += 1
                continue
            if "PetriSpot" in line and "Exec failed" in line:
                petrispot_fail += 1
                continue
            if line.startswith("PetriSpot timed out after"):
                petrispot_timeout += 1
                continue
            m = RE_WALK.match(line)
            if m:
                walk[0] += 1
                walk[1] += int(m.group(2))
                walk[2] += int(m.group(1))
                walk[3] += int(m.group(3))
                walk[4] += int(m.group(4)) != 0
                continue
            if not failure:
                for name, pattern in FAILURES:
                    if pattern in line:
                        failure = name
                        break
            if not run["start"]:
                m = RE_STAMP.match(line)
                if m:
                    run["start"] = m.group(1)
            m = RE_TC.search(line)
            if m:
                kind, name, duration = m.group(1), m.group(2), m.group(3)
                if kind == "testStarted":
                    started += 1
                elif kind == "testFailed":
                    failed += 1
                elif kind == "testFinished":
                    finished += 1
                    if name == "all":
                        run["duration(ms)"] = duration or ""
                    elif name == "runits":
                        run["first(ms)"] = duration or ""
                continue
            m = RE_HOST.match(line)
            if m:
                run["host"] = m.group(1)
                continue
            m = RE_TIMEOUT.match(line)
            if m:
                run["timeout(s)"] = m.group(1)
                continue
            m = RE_VERSION.match(line)
            if m:
                run["version"] = m.group(1)
                continue
            m = RE_TITLE.match(line)
            if m and not run["Model"]:
                run["Examination"] = m.group(1).rsplit("-", 1)[-1]

    verdicts = []
    ok = wrong = missed = bonus = extra = known = 0
    for name in sorted(set(oracle) | set(answer)):
        exp = oracle.get(name)
        got = answer.get(name)
        if exp is not None and exp != "?":
            known += 1
        if exp is None:
            status, extra = "extra", extra + 1
        elif got is None:
            if exp == "?":
                status = "none"
            else:
                status, missed = "missed", missed + 1
        elif exp == "?":
            status, bonus = "bonus", bonus + 1
        elif equivalent(exp, got):
            status, ok = "ok", ok + 1
        else:
            status, wrong = "wrong", wrong + 1
        verdicts.append({
            "log": path,
            "Model": run["Model"],
            "Examination": run["Examination"],
            "formula": name,
            "oracle": "" if exp is None else exp,
            "answer": "" if got is None else got,
            "status": status,
            "techniques": techniques.get(name, ""),
        })

    run.update({
        "petrispot": petrispot, "petrispot fail": petrispot_fail,
        "petrispot timeout": petrispot_timeout,
        "walk calls": walk[0], "walk asked": walk[1], "walk solved": walk[2],
        "walk ms": walk[3], "walk killed": walk[4], "failure": failure, "Test started": started, "Test fail": failed, "Test fin": finished,
        "oracle": len(oracle), "known": known, "answered": len(answer),
        "ok": ok, "wrong": wrong, "missed": missed, "bonus": bonus,
        "extra": extra,
    })
    return run, verdicts


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("dirs", nargs="+", help="directories holding *out logs")
    ap.add_argument("-o", "--outdir", default=".", help="where to write the CSV")
    ap.add_argument("--prefix", default="", help="prefix of the CSV names")
    args = ap.parse_args()

    logs = []
    for d in args.dirs:
        if os.path.isdir(d):
            logs += sorted(os.path.join(d, f) for f in os.listdir(d)
                           if f.endswith("out"))
        else:
            logs.append(d)

    os.makedirs(args.outdir, exist_ok=True)
    runs_path = os.path.join(args.outdir, args.prefix + "runs.csv")
    verdicts_path = os.path.join(args.outdir, args.prefix + "verdicts.csv")
    with open(runs_path, "w", newline="") as rf, \
         open(verdicts_path, "w", newline="") as vf:
        rw = csv.DictWriter(rf, RUN_COLUMNS)
        vw = csv.DictWriter(vf, VERDICT_COLUMNS)
        rw.writeheader()
        vw.writeheader()
        for path in logs:
            run, verdicts = parse(path)
            rw.writerow(run)
            vw.writerows(verdicts)
    print(f"{len(logs)} logs -> {runs_path}, {verdicts_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
