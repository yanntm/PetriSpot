#!/usr/bin/env python3
"""Does the per property reduction phase of UpperBounds pay for itself?

`Application.java` runs a global `UpperBoundsSolver.applyReductions` over the
whole cohort, and when values remain it loops once per unsolved property,
copying the net, keeping that property alone and reducing again -- the phase
the log brackets with

    Starting property specific reduction for <formula>
    Ending property specific reduction for <formula> in <t> ms.

For every such bracket this counts the seconds spent and the FORMULA verdicts
that appeared inside it. Run from a directory of examination folders:

    perprop_payoff.py UB
"""

import os
import re
import sys
from collections import Counter

RE_START = re.compile(r"^Starting property specific reduction for (\S+)")
RE_END = re.compile(r"^Ending property specific reduction for (\S+) in (\d+) ms\.")
RE_FORMULA = re.compile(r"^FORMULA (\S+) (\S+)")


def main(dirs):
    logs = reached = 0
    brackets = solved_inside = solved_self = 0
    ms_total = ms_paying = 0
    per_model = Counter()
    for d in dirs:
        for name in sorted(os.listdir(d)):
            if not name.endswith("out"):
                continue
            logs += 1
            path = os.path.join(d, name)
            current = None
            found = []
            seen_phase = False
            for line in open(path, errors="replace"):
                m = RE_START.match(line)
                if m:
                    seen_phase = True
                    current, found = m.group(1), []
                    continue
                m = RE_END.match(line)
                if m and current is not None:
                    ms = int(m.group(2))
                    brackets += 1
                    ms_total += ms
                    if found:
                        solved_inside += 1
                        ms_paying += ms
                        per_model[os.path.basename(d) + " " + m.group(1).rsplit("-", 2)[0]] += 1
                        if current in found:
                            solved_self += 1
                    current = None
                    continue
                if current is not None:
                    m = RE_FORMULA.match(line)
                    if m:
                        found.append(m.group(1))
            reached += seen_phase

    print(f"logs                     : {logs}")
    print(f"reached the phase        : {reached}")
    print(f"per property reductions  : {brackets}")
    print(f"  ... that produced any verdict : {solved_inside}")
    print(f"  ... for the very property they targeted : {solved_self}")
    print(f"time in the phase        : {ms_total/1000:.0f} s, of which {ms_paying/1000:.0f} s produced something")
    if brackets:
        print(f"payoff rate              : {100*solved_inside/brackets:.1f} %")
    for k, v in per_model.most_common(10):
        print(f"    {v:3d}  {k}")


if __name__ == "__main__":
    main(sys.argv[1:])
