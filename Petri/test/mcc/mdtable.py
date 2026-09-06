"""Markdown helpers shared by the report writers: tables, counts, families."""

import re

RE_FAMILY = re.compile(r"-(PT|COL)-.*$")


def I(x):
    """An integer column, empty cells being 0."""
    return int(x or 0)


def family(model):
    """`ResIsolation` for `ResIsolation-PT-N12P1`."""
    return RE_FAMILY.sub("", model)


def table(headers, rows, align=None):
    """A markdown table; align is a string of 'l'/'r' per column (right for numbers by default)."""
    if align is None:
        align = "".join("r" if rows and isinstance(rows[0][i], (int, float)) else "l" for i in range(len(headers)))
    out = ["| " + " | ".join(str(h) for h in headers) + " |",
           "|" + "|".join(" ---: " if a == "r" else " --- " for a in align) + "|"]
    for r in rows:
        out.append("| " + " | ".join(fmt(c) for c in r) + " |")
    return "\n".join(out) + "\n"


def fmt(c):
    if isinstance(c, float):
        return f"{c:.2f}" if abs(c) < 100 else f"{c:.0f}"
    return str(c)


def hours(ms):
    return I(ms) / 3.6e6


def seconds(ms):
    return I(ms) / 1000


def at_wall(run, timeout_s):
    """A run that was killed by the harness, with a little slack for the trailer."""
    return I(run.get("duration(ms)")) > (timeout_s - 50) * 1000


def truncated(run):
    """A log without its trailer: no duration at all (killed by the scheduler, still running)."""
    return not I(run.get("duration(ms)"))


def bullet(items, cap=40):
    """A bulleted list, elided past cap items."""
    out = "".join(f"* {it}\n" for it in items[:cap])
    if len(items) > cap:
        out += f"* ... {len(items) - cap} more\n"
    return out
