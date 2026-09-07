#!/usr/bin/env python3
"""The Liveness examination of a PNML as CTL: one `(ctl liveI (AG (EF (fireable tI))))` property per
transition, or with --one the single conjunction `(ctl Liveness (and (AG (EF (fireable t)))...))`, which
is FALSE as soon as one transition is not live.

    liveness_props.py [--one] model.pnml out.sexpr [N]      # the first N transitions in definition order (default all)
"""

import re
import sys

args = sys.argv[1:]
one = "--one" in args
args = [a for a in args if a != "--one"]
pnml, out = args[0], args[1]
n = int(args[2]) if len(args) > 2 else None
names = re.findall(rb'<transition\s[^>]*?id="([^"]+)"', open(pnml, "rb").read())
if n is not None:
    names = names[:n]
with open(out, "w") as f:
    if one:
        f.write("(ctl Liveness (and\n")
        for t in names:
            f.write(f"  (AG (EF (fireable {t.decode()})))\n")
        f.write("))\n")
    else:
        for i, t in enumerate(names):
            f.write(f"(ctl live{i} (AG (EF (fireable {t.decode()}))))\n")
print(f"{len(names)} transitions -> {out}", file=sys.stderr)
