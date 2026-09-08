# LTSmin partial order reduction misses an accepting cycle

`pins2lts-seq` under `-p` reports the product empty on a model where the same
binary without `-p` finds an accepting cycle. The verdict ITS-Tools publishes
is then TRUE where the property is FALSE.

Seen on `StigmergyCommit-PT-02b-LTLCardinality-03` in campaign 2026-09-07 (one
wrong verdict in 30 373), and reproduced with the products `202609071637` and
`202609080210` alike -- it is not a regression, it is a race: whichever engine
answers first wins, and the stuttering random walk usually finds the witness
before LTSmin is asked. When it does not, the wrong verdict is published.

The formula is `!(G((F(p0)||X(X(F(p1))))))`, the negation ITS-Tools checks for
emptiness; the automaton is very weak, state based and stutter invariant, so
partial order reduction is applied.

## Reproduce

`model.c` is the PINS model ITS-Tools generates for the instance (43 places,
107 transitions after the SI_LTL reductions), `stateBased.hoa` the automaton
handed to LTSmin, `curaut.hoa` the automaton it was derived from
(`autfilt --buchi --state-based-acceptance --deterministic curaut.hoa`).

```
B=<product>/plugins/fr.lip6.move.gal.ltsmin.binaries_*/bin
gcc -c -I$B/include/ -I. -std=c99 -fPIC -O0 model.c && gcc -shared -o gal.so model.o
$B/pins2lts-seq-linux64 ./gal.so -p --pins-guards --when --hoa stateBased.hoa --buchi-type=spotba
$B/pins2lts-seq-linux64 ./gal.so             --when --hoa stateBased.hoa --buchi-type=spotba
```

```
with -p     : Empty product with LTL!      state space 101 levels, 8789 states 25220 transitions
without -p  : Accepting cycle FOUND!
```

`--proviso=color` and `--proviso=stack` (the default here) both report empty,
and `--pins-guards` changes nothing.

## Visibility is correct; the reduction loses the cycle anyway

Every reduced run prints `Visible groups: 0 / 111, labels: 2 / 111`, which is a
trap: `pins2lts-seq.c` prints the `GBgetPorGroupVisibility` **array**, which only
`pins_add_group_visible` ever writes. The reduction works from `ctx->visible`, a
separate `bms_t` that `init_visible_labels` (`pins2pins-por.c`) fills, per state,
by unioning the NES and NDS groups of every visible label.

Instrumenting that function settles it (upstream at 07f9bf8, built with Spot):

```
Initializing POR dependencies: labels 109, guards 107
Visible groups: 0 / 111, labels: 2 / 111
PORDIAG before: 2 visible labels,  0 visible groups
PORDIAG after:  2 visible labels, 22 visible groups
Empty product with LTL!
```

POR receives exactly the 22 groups our model declares, and still reports the
product empty where the same binary without `-p` finds an accepting cycle. The
defect is in the reduction, downstream of visibility.

Everything the model owes LTSmin is there:

* `label_matrix` gives label 107 (`LTLAPp0`) the state variables 12 and 28, and
  label 108 (`LTLAPp1`) the variable 25.
* `write_matrix` has 22 groups writing one of those three:
  `[1, 2, 3, 6, 7, 28, 29, 30, 31, 32, 40, 46, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63]`.
* `mayEnableAtom` and `mayDisableAtom`, the NES and NDS rows of those two labels,
  list exactly those groups (20 for `p0`, 2 for `p1`) -- and those are what
  `init_visible_labels` turns into the 22 above.

Two further candidates are ruled out. We never emit `GBsetDMInfoMayWrite`
(`Gal2PinsTransformerNext`, `PetriNet2PinsTransformer`) though POR's dependency
pruning reads may-write, but adding it changes nothing and LTSmin defaults
may-write to the combined matrix anyway (`pins.c`), a superset. And the property
is genuinely stutter invariant, so `-p` is legitimately applicable: Spot on the
source formula gives `properties: stutter-invariant very-weak weak
inherently-weak`.

What remains is an asymmetry in LTSmin itself. Its LTL front end refuses the
combination outright,

```
pins2lts-seq ./gal.so -p --pins-guards --when \
   --ltl='[]((<>(LTLAPp0==true)) || (X (X (<>(LTLAPp1==true)))))' --ltl-semantics=spin
** error **: The neXt operator is not allowed in combination with --por
```

while the same front end without `-p` finds the accepting cycle. That guard
(`pins2pins-ltl.c`) is syntactic -- it refuses any X, though this formula is
stutter invariant -- and it lives only on the path where LTSmin parses a formula.
Through `--hoa` it sees an automaton, the guard cannot fire, and the reduction
runs. `-m` shows POR being set up before the automaton exists at all:
`Initializing POR dependencies: labels 109, guards 107`, our raw counts, printed
before `buchi has 3 states`.

Reproducing the instrumentation: `~/git/ltsmin` is upstream
`utwente-fmt/ltsmin`; `--hoa` needs Spot at configure time
(`PKG_CONFIG_PATH=/usr/local/lib/pkgconfig`, `SPOT = yes`), the bundled `lemon`
needs `CFLAGS=-std=gnu17` under GCC 15, and because the build is configured
`--disable-dependency-tracking` a `make clean` is required after re-running
configure or the objects keep the old `config.h`.

## Every reduction variant loses it, and the witness is shallow

`--por` takes an algorithm and the search a proviso. Both working algorithms, on
all three provisos, report the product empty:

| `--por` | `--proviso` | verdict | reduced state space |
| --- | --- | --- | --- |
| heur | stack | Empty product | 8789 states, 25220 transitions |
| heur | color | Empty product | 8789 states, 25220 transitions |
| heur | closedset | Empty product | 8533 states, 24517 transitions |
| del | stack | Empty product | 9125 states, 26045 transitions |
| del | color | Empty product | 9125 states, 26045 transitions |
| del | closedset | Empty product | 8869 states, 25342 transitions |

(`tr` and `str` abort with "Undefined PC identification criteria", they want a
transaction structure this model has none of.) So it is not a proviso quirk and
not one stubborn set heuristic: every reduction here drops the cycle, while the
unreduced search finds it in 5 ms, inside an SCC at depth 1039, and writes a
70 step witness:

```
pins2lts-seq ./gal.so --when --hoa stateBased.hoa --buchi-type=spotba --trace=cex.gcf
ltsmin-printtrace cex.gcf          # length of trace is 70
```

The automaton is weak, and the LTL layer says so:
`Weak Buchi automaton detected, adding non-accepting as progress label`
(`pins2pins-ltl.c`, `is_weak(ba)` for BA and SPOTBA types, adding
`LTSMIN_STATE_LABEL_WEAK_LTL_PROGRESS` as label 110 beside accepting at 109).
Our own labels stay at 107 and 108. That weak-automaton progress path is the
next thing to look at, being the one piece of machinery specific to this shape
of automaton.

## It is the V proviso, and `--no-V` avoids it

POR keeps hidden switches for its provisos. Swapping the visibility proviso
restores the right answer with the reduction still on:

| `--por` | flags | verdict |
| --- | --- | --- |
| heur | (default) | Empty product -- wrong |
| heur | `--no-L12`, `--no-mc`, `--no-mcnds` | Empty product -- wrong |
| heur | `--no-V` | Accepting cycle FOUND |
| heur | `--weak` | Accepting cycle FOUND |
| del | any of the above | Empty product -- wrong |

`--no-V` replaces LTSmin's own visibility proviso with Peled's; `--weak` swaps
the stubborn set theory. Both make this case sound, which puts the defect in the
default V proviso. The deletion algorithm `--por=del` is wrong under every
combination and should not be used at all.

What it costs, measured on a product that really is empty so the whole space is
walked (`--hoa curaut.hoa --buchi-type=tgba`):

| setting | states | transitions |
| --- | ---: | ---: |
| default | 8789 | 25220 |
| `--no-V` | 10099 | 30013 |
| `--weak` | 12146 | 38176 |
| no POR | 12183 | 39656 |

So `--no-V` keeps most of the reduction, a sixth off the unreduced space, while
`--weak` gives up nearly all of it. `LTSminRunner.checkProperty` already carries
the line commented out beside the `-p` it passes:

```java
ltsmin.addArg("-p");
ltsmin.addArg("--pins-guards");
//ltsmin.addArg("--no-V");
```

Uncommenting it is the fix on our side, and needs no patched LTSmin. It is not a
proof of soundness everywhere -- Peled's proviso is the textbook one and LTSmin's
optimised replacement is what misbehaves here -- but it turns a wrong answer into
a correct one at a measured cost.

Ruled out along the way, none of them the cause: the weak Buchi progress label
(`ctx->is_weak` forced false changes nothing), the `SAFETY` flag (set after
visibility, and every use is `SAFETY || PINS_LTL` with PINS_LTL true), the cycle
proviso choice, `--no-mc`, `--no-mcnds`, and the may-write matrix.

## What to do with it

Ours to work around, upstream to fix. Our side of the choice is
`LTSminRunner.checkProperty`, where `if (doPOR && isStutterInvariant(pbody))`
adds `-p`: the condition is right in theory and the reduction is wrong in
practice, so the workarounds are to withhold `-p` when the formula carries an X
(what LTSmin does for itself on the `--ltl` path), or to trust a witness over an
emptiness claim when the two engines disagree. The reproduction folder built from the MCC archive is
`bench/models/Stigmergy-f03/` (git ignored), made by
`Petri/test/extract_property.py`; `Petri/test/logs/stig-dbg-1.log` is a full
run that takes the wrong branch, `Petri/test/logs/stigmergy-f03-debug.log` one
that answers FALSE by `STUTTER_TEST` before LTSmin is reached.
