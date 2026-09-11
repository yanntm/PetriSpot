#!/usr/bin/env python3
"""Check one free-cycle net through native and standalone counting reductions."""
from __future__ import annotations

import argparse
import math
from pathlib import Path
import subprocess
import struct
import tempfile
import time


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hsc", type=Path, required=True)
    parser.add_argument("--petri", type=Path, required=True)
    parser.add_argument("--places", type=int, default=2)
    parser.add_argument("--tokens", type=int, default=3)
    parser.add_argument("--components", type=int, default=1)
    parser.add_argument("--skip-original", action="store_true",
                        help="use only the analytic oracle for large Cartesian products")
    args = parser.parse_args()
    if args.places < 2 or args.tokens < 0 or args.components < 1:
        parser.error("need >=2 places, >=0 tokens and >=1 components")
    logs = Path(__file__).resolve().parents[1] / "logs"
    logs.mkdir(exist_ok=True)
    tag = f"counting-chain-{args.places}-{args.tokens}-{args.components}"
    expected = {
        "STATES": str(math.comb(args.tokens + args.places - 1, args.places - 1) ** args.components),
        "MAX_TOKEN_IN_PLACE": str(max(2, args.tokens)),
        "MAX_TOKEN_PER_MARKING": str(2 + args.tokens * args.components),
    }
    deadline = time.monotonic() + 15
    with (logs / f"{tag}.log").open("w") as log, tempfile.TemporaryDirectory(dir=logs) as tmp:
        root = Path(tmp)
        model = root / "model.pnml"
        forms = ['<pnml><net id="cycle" type="http://www.pnml.org/version-2009/grammar/ptnet"><page id="page">',
                 '<place id="constant"><initialMarking><text>2</text></initialMarking></place>']
        for c in range(args.components):
            for p in range(args.places):
                forms.append(f'<place id="p{c}_{p}"><initialMarking><text>{args.tokens if p == 0 else 0}</text></initialMarking></place>')
                forms.append(f'<transition id="t{c}_{p}"/>')
                forms.append(f'<arc id="a{c}_{p}" source="p{c}_{p}" target="t{c}_{p}"/>')
                forms.append(f'<arc id="b{c}_{p}" source="t{c}_{p}" target="p{c}_{(p + 1) % args.places}"/>')
        model.write_text("\n".join(forms + ["</page></net></pnml>"]))

        def run(command: list[str], label: str, check: bool = True) -> None:
            result = subprocess.run(command, capture_output=True, text=True,
                                    timeout=max(0.01, deadline - time.monotonic()))
            log.write(f"== {label}\n{result.stdout}{result.stderr}\n")
            log.flush()
            if result.returncode:
                raise AssertionError(f"{label}: exit {result.returncode}; see {log.name}")
            if check:
                values = {words[1]: words[2] for line in result.stdout.splitlines()
                          if (words := line.split()) and words[0] == "STATE_SPACE" and len(words) >= 3}
                for key, value in expected.items():
                    if values.get(key) != value:
                        raise AssertionError(f"{label}: {key} expected {value}, got {values.get(key)}; see {log.name}")
                counts = [line for line in result.stderr.splitlines() if line.startswith("hsc-pn: counts ")]
                if counts:
                    fields = dict(word.split("=", 1) for word in counts[-1].split()[2:])
                    assert fields["reach_weighted_states"] == expected["STATES"], fields
                print(f"{label}: {values}")

        def check_export(path: Path) -> None:
            data = path.read_bytes()
            assert struct.unpack_from("<I", data, 6)[0] == 0, "constant components remain in net"
            if args.tokens > 0:
                assert b"PCONST\0\0" in data, "missing constant component record"

        hsc = [str(args.hsc.resolve()), "--states", "--totalTime", "3", "-v"]
        if not args.skip_original:
            run(hsc + ["-i", str(model)], "original")
        native = root / "native.pnet"
        run(hsc + ["-i", str(model), "--reduce", "--export-net", str(native)], "native reduce")
        check_export(native)
        run(hsc + ["--net", str(native), "--reduce"], "native re-reduce")
        previous = model
        for iteration in range(2):
            target = root / f"export-{iteration}"
            run([str(args.petri.resolve()), "reduce", "-i" if iteration == 0 else "--net",
                 str(previous), "--goal", "STATESPACE", "--deadMs", "0", "--output", str(target)],
                f"standalone export {iteration}", check=False)
            previous = target / "model.pnet"
            check_export(previous)
            run(hsc + ["--net", str(previous)], f"standalone count {iteration}")
    print(f"PASS {tag}")


if __name__ == "__main__":
    main()
