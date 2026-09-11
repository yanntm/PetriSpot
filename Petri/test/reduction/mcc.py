#!/usr/bin/env python3
"""Local real-CLI reduction validation; one model subprocess, hard 15 s total."""
from __future__ import annotations

import argparse
import fnmatch
import json
from pathlib import Path
import re
import subprocess
import sys
import tarfile
import tempfile
import time
from typing import Any

EXAMS = {"RC": "ReachabilityCardinality", "RF": "ReachabilityFireability",
         "UB": "UpperBounds", "CTLC": "CTLCardinality", "CTLF": "CTLFireability",
         "RD": "ReachabilityDeadlock"}


def verdicts(text: str) -> dict[str, str]:
    return dict(re.findall(r"^FORMULA\s+(\S+)\s+(\S+)", text, re.M))


def emit(row: dict[str, Any]) -> None:
    print(json.dumps(row, sort_keys=True), flush=True)


def one_model(args: argparse.Namespace) -> None:
    archive = Path(args.one)
    model = archive.name.removesuffix(".tgz")
    root = Path(args.workdir)
    # Extract only the input and query files; never write archive links/paths.
    wanted = {"model.pnml", *(name + ".xml" for name in EXAMS.values())}
    with tarfile.open(archive) as source:
        for member in source:
            leaf = Path(member.name).name
            if member.isfile() and leaf in wanted:
                stream = source.extractfile(member)
                assert stream is not None
                with (root / leaf).open("wb") as output:
                    while block := stream.read(1024 * 1024):
                        output.write(block)
    for exam in args.exams.split(","):
        property_file = root / (EXAMS[exam] + ".xml")
        if exam != "RD" and not property_file.exists():
            emit({"model": model, "exam": exam, "status": "missing-properties"})
            continue
        oracle_file = Path(args.oracles) / f"{model}-{exam}.out"
        oracle = verdicts(oracle_file.read_text()) if oracle_file.exists() else {}
        pair: dict[str, dict[str, str]] = {}
        for mode in ("reduced", "original"):
            command = [str(Path(args.binary).resolve()), "-i", str(root / "model.pnml"),
                       "--totalTime=1", "-t", "1", "--sweepTime=0", "--walkSteps=100",
                       "--tasks=1", "--threads=1", "--spawn=off", "--slice=128",
                       "--runLength=100", "--strategy=random", "--seed=1", "--printUnknown",
                       "--ctlRounds=1", "--ctlSteps=100", "--ctlRegion=100", "--ctlRunLength=100", "-q"]
            if exam == "RD":
                command.append("--findDeadlock")
            else:
                command.append("--props=" + str(property_file))
            reduction_text = ""
            if mode == "reduced" and args.standalone:
                destination = root / ("reduced-" + exam)
                if exam == "RD":
                    property_file = root / "deadlock.sexpr"
                    property_file.write_text('(deadlock ReachabilityDeadlock)\n')
                remaining = 14.5 - (time.monotonic() - args.started)
                if remaining < 0.1:
                    emit({"model": model, "status": "model-timeout"})
                    return
                transform = subprocess.run([
                    "timeout", str(remaining) + "s", str(Path(args.binary).resolve()), "reduce",
                    "-i", str(root / "model.pnml"), "--props", str(property_file),
                    "--output", str(destination), "--reductionMs=10000"],
                    capture_output=True, text=True, check=False)
                reduction_text = transform.stdout + transform.stderr
                if transform.returncode:
                    emit({"model": model, "exam": exam, "mode": mode,
                          "status": "timeout" if transform.returncode == 124 else "error",
                          "error_tail": reduction_text[-4000:]})
                    return
                command[1:3] = ["--net", str(destination / "model.pnet")]
                command = [arg for arg in command if not arg.startswith("--props=") and arg != "--findDeadlock"]
                command.append("--props=" + str(destination / "properties.sexpr"))
            elif mode == "reduced":
                command += ["--reduce", "--reductionMs=10000"]
            if args.lp and exam in ("RC", "RF", "UB"):
                command += ["--lp", "--lpTime=0.02", "--lpSolves=50"]
            start = time.monotonic()
            # Inner timeout guarantees no orphaned model process if the outer
            # model worker reaches its own limit. GNU timeout remains its parent.
            remaining = 14.5 - (time.monotonic() - args.started)
            if remaining < 0.1:
                emit({"model": model, "status": "model-timeout"})
                return
            with tempfile.TemporaryFile(mode="w+", dir=args.logs) as log:
                completed = subprocess.run(["timeout", str(remaining) + "s", *command],
                                           stdout=log, stderr=subprocess.STDOUT, check=False)
                log.seek(0)
                text = reduction_text + log.read()
            answers = verdicts(text)
            if exam == "RD" and "ReachabilityDeadlock" in answers and len(oracle) == 1:
                answers[next(iter(oracle))] = answers.pop("ReachabilityDeadlock")
            wrong = {name: {"actual": answer, "oracle": oracle[name]}
                     for name, answer in answers.items()
                     if name in oracle and oracle[name] != "?" and oracle[name] != answer}
            bounds = dict(re.findall(r"^BOUND\s+(\S+)\s+(-?\d+)", text, re.M))
            for name, bound in bounds.items():
                expected = oracle.get(name, "?")
                if expected.isdigit() and int(bound) > int(expected):
                    wrong[name] = {"lower_bound": bound, "oracle": expected}
            pair[mode] = answers
            row: dict[str, Any] = {
                "model": model, "exam": exam, "mode": mode, "returncode": completed.returncode,
                "seconds": round(time.monotonic() - start, 4), "answers": answers,
                "oracle_matches": sum(oracle.get(n) == v for n, v in answers.items()),
                "unverified_answers": sum(oracle.get(n, "?") == "?" for n in answers),
                "wrong": wrong, "bounds": bounds,
                "reduction": re.findall(r"^Reduction .*", text, re.M),
                "status": "ok" if completed.returncode == 0 else "timeout" if completed.returncode == 124 else "error"}
            if completed.returncode != 0:
                row["error_tail"] = text[-4000:]
            if mode == "original":
                row["conflicts"] = {n: [v, pair["reduced"][n]] for n, v in answers.items()
                                    if n in pair["reduced"] and pair["reduced"][n] != v}
            emit(row)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", default="build/petri64")
    parser.add_argument("--inputs", required=True)
    parser.add_argument("--oracles", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--model", default="*-PT-*")
    parser.add_argument("--max-models", type=int, default=0)
    parser.add_argument("--exams", default=",".join(EXAMS))
    parser.add_argument("--lp", action="store_true")
    parser.add_argument("--standalone", action="store_true", help="Export a reduced pair, then solve it in a separate process")
    parser.add_argument("--one", help=argparse.SUPPRESS)
    parser.add_argument("--workdir", help=argparse.SUPPRESS)
    parser.add_argument("--logs", default="Petri/test/logs")
    args = parser.parse_args()
    args.started = time.monotonic()
    if args.one:
        one_model(args)
        return 0
    archives = sorted(p for p in Path(args.inputs).glob("*.tgz")
                      if fnmatch.fnmatch(p.stem, args.model))
    if args.max_models:
        archives = archives[:args.max_models]
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    Path(args.logs).mkdir(parents=True, exist_ok=True)
    totals: dict[str, int] = {"models": 0, "runs": 0, "wrong": 0, "conflicts": 0,
                              "errors": 0, "timeouts": 0, "oracle_matches": 0}
    with output.open("x") as results:
        for archive in archives:
            # Parent owns cleanup even when the worker reaches its hard deadline.
            with tempfile.TemporaryDirectory(prefix=".reduction-", dir=args.inputs) as directory:
                command = [sys.executable, __file__, *sys.argv[1:], "--one", str(archive),
                           "--workdir", directory]
                try:
                    run = subprocess.run(command, capture_output=True, text=True, timeout=15, check=False)
                    text = run.stdout
                    if run.returncode:
                        text += json.dumps({"model": archive.stem, "status": "error",
                                            "error_tail": run.stderr[-4000:]}) + "\n"
                except subprocess.TimeoutExpired as error:
                    text = error.stdout or b""
                    if isinstance(text, bytes):
                        text = text.decode(errors="replace")
                    text += json.dumps({"model": archive.stem, "status": "model-timeout"}) + "\n"
                for line in text.splitlines():
                    row = json.loads(line)
                    results.write(json.dumps(row, sort_keys=True) + "\n")
                    totals["runs"] += "mode" in row
                    totals["wrong"] += len(row.get("wrong", {}))
                    totals["conflicts"] += len(row.get("conflicts", {}))
                    totals["errors"] += row.get("status") == "error"
                    totals["timeouts"] += row.get("status") in ("timeout", "model-timeout")
                    totals["oracle_matches"] += row.get("oracle_matches", 0)
            totals["models"] += 1
            results.flush()
            if totals["models"] % 25 == 0:
                print(json.dumps(totals), flush=True)
    print(json.dumps(totals), flush=True)
    return int(totals["wrong"] > 0 or totals["conflicts"] > 0 or totals["errors"] > 0)


if __name__ == "__main__":
    raise SystemExit(main())
