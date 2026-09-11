#!/usr/bin/env python3
"""Summarize a native reduction validation JSONL without treating UNKNOWN as agreement."""
from __future__ import annotations
import collections
import json
from pathlib import Path
import re
import sys


def main() -> None:
    counts: collections.Counter[str] = collections.Counter()
    rules: collections.Counter[str] = collections.Counter()
    models: set[str] = set()
    for line in Path(sys.argv[1]).read_text().splitlines():
        row = json.loads(line)
        models.add(row["model"])
        mode = row.get("mode", "worker")
        counts[f"{mode}.{row['status']}"] += 1
        counts[f"{mode}.answers"] += len(row.get("answers", {}))
        counts[f"{mode}.oracle_matches"] += row.get("oracle_matches", 0)
        counts[f"{mode}.unverified"] += row.get("unverified_answers", 0)
        counts[f"{mode}.wrong"] += len(row.get("wrong", {}))
        counts["conflicts"] += len(row.get("conflicts", {}))
        for reduction in row.get("reduction", []):
            match = re.match(r"Reduction rule (.*?): (\d+) sparse edits", reduction)
            if match:
                rules[match[1]] += int(match[2])
            if "(limit)" in reduction:
                counts["reduction.limit"] += 1
        if row.get("wrong") or row.get("conflicts") or row["status"] == "error":
            print(json.dumps(row, sort_keys=True))
    print(json.dumps({"models": len(models), "counts": dict(counts), "rule_edits": dict(rules)}, sort_keys=True))


if __name__ == "__main__":
    main()
