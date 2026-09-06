#!/usr/bin/env python3
"""One `(reach pI (fireable tI))` target per transition of a PNML, the QuasiLivenessAll query set.

    qla_props.py model.pnml out.sexpr [N]      # the first N transitions in definition order (default all)

Mirrors what ITS-Tools hands PetriSpot for QuasiLivenessAll, to measure the
walker alone on a large target set.
"""

import re
import sys

pnml, out = sys.argv[1], sys.argv[2]
n = int(sys.argv[3]) if len(sys.argv) > 3 else None
names = re.findall(rb'<transition\s[^>]*?id="([^"]+)"', open(pnml, "rb").read())
if n is not None:
    names = names[:n]
with open(out, "w") as f:
    for i, t in enumerate(names):
        f.write(f"(reach prop{i} (fireable {t.decode()}))\n")
print(f"{len(names)} targets -> {out}", file=sys.stderr)
