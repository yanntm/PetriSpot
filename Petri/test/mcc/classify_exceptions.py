#!/usr/bin/env python3
"""classify_exceptions.py: what went wrong in a collected examination folder.

    python3 Petri/test/mcc/classify_exceptions.py /data/ythierry/MCC26run/<date>/<EXAM> [-v]

One pass over the OAR.*.stdout logs. Each log is attributed the exception families it
carries (the message shape, numbers folded), its stack trace size in frames, and how many
FORMULA verdicts it produced. Prints the families by log count, then the logs that print a
stack trace, largest first: a raw trace is a failure nobody caught and named, and a native
closed-world miss (MissingReflectionRegistrationError) hides there.
"""
import re, sys, os, glob, collections

FAMILY = [
    (r"MissingReflectionRegistrationError", "native closed world: MissingReflectionRegistrationError"),
    (r"Could not find partition element", "partition refinement: could not find partition element"),
    (r"OutOfMemoryError", "OutOfMemoryError"),
    (r"OverlargeMarkingException", "OverlargeMarkingException (annotations too wide)"),
    (r"TimeoutException", "TimeoutException on a subprocess"),
    (r"IOException: (Broken pipe|Stream closed)", "IOException broken pipe / stream closed"),
    (r'Timeout of Z3 solver reached', "Z3 solver timeout"),
    (r"NullPointerException", "NullPointerException"),
    (r"IndexOutOfBoundsException", "IndexOutOfBoundsException"),
    (r"([A-Za-z0-9_.]*(?:Exception|Error))\b", None),   # catch-all, keeps the class name
]
AT = re.compile(r"^\s*at \S+\(")
TEST = re.compile(r"Running test : oracle\.(\S+)")

def families(text):
    out = set()
    for line in text.splitlines():
        if "Exception" not in line and "Error" not in line:
            continue
        for pat, name in FAMILY:
            m = re.search(pat, line)
            if m:
                out.add(name or m.group(1))
                break
    return out

def main(folder, verbose=False):
    rows = []
    for f in sorted(glob.glob(os.path.join(folder, "*.stdout"))):
        text = open(f, errors="replace").read()
        m = TEST.search(text)
        rows.append(dict(log=os.path.basename(f), test=m.group(1) if m else "?",
                         fam=families(text),
                         frames=sum(1 for l in text.splitlines() if AT.match(l)),
                         verdicts=text.count("\nFORMULA ") + text.startswith("FORMULA ")))
    print(f"{len(rows)} logs, {sum(r['verdicts'] for r in rows)} verdicts, "
          f"{sum(1 for r in rows if r['verdicts'] == 0)} silent, "
          f"{sum(1 for r in rows if r['frames'])} printing a stack trace\n")
    count = collections.Counter(fam for r in rows for fam in r["fam"])
    print("family                                                    logs  silent  verdicts")
    for fam, n in count.most_common():
        hit = [r for r in rows if fam in r["fam"]]
        print(f"{fam:55.55s} {n:5d} {sum(1 for r in hit if r['verdicts']==0):7d} {sum(r['verdicts'] for r in hit):9d}")
    print("\nlogs printing a stack trace, largest first")
    print("frames  verdicts  test                                        families")
    for r in sorted((r for r in rows if r["frames"]), key=lambda r: -r["frames"])[:40 if not verbose else 10**9]:
        print(f"{r['frames']:6d} {r['verdicts']:9d}  {r['test']:43.43s} {', '.join(sorted(r['fam']))}")

if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "-v"]
    main(args[0], "-v" in sys.argv)
