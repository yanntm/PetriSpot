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

## What we know, and what the `Visible groups: 0` line is not

Every reduced run prints

```
Visible groups: 0 / 111, labels: 2 / 111
```

That line is the obvious suspect and it is a trap. `pins2lts-seq.c` prints the
`GBgetPorGroupVisibility` **array**, which only `pins_add_group_visible` ever
writes. The POR layer derives its working set elsewhere: `init_visible_labels`
(`pins2pins-por.c`) unions, for every visible label, the groups of that label's
NES and NDS into `ctx->visible`, a separate `bms_t`. It is called from
`por_init_transitions`, per state during exploration, so it runs long after the
HOA layer has parsed the atomic propositions. A zero in the printed array is
therefore consistent with POR having the right groups internally. We have not
shown that visibility is empty where it matters.

What is established is that the model we hand over carries the information:

* `label_matrix` gives label 107 (`LTLAPp0`) the state variables 12 and 28, and
  label 108 (`LTLAPp1`) the variable 25.
* `write_matrix` has 22 groups writing one of those three:
  `[1, 2, 3, 6, 7, 28, 29, 30, 31, 32, 40, 46, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63]`.
* `mayEnableAtom` and `mayDisableAtom`, the NES and NDS rows for those two
  labels, list exactly those groups (20 for `p0`, 2 for `p1`) -- precisely what
  `init_visible_labels` consumes.

And two candidate faults on our side are ruled out:

* We register `GBsetDMInfoRead`, `GBsetDMInfoMustWrite` and `GBsetDMInfo`, never
  `GBsetDMInfoMayWrite` (`Gal2PinsTransformerNext`, `PetriNet2PinsTransformer`),
  and POR's dependency pruning reads may-write. Adding
  `GBsetDMInfoMayWrite(m, wm)` to `model.c` and rebuilding changes nothing.
  LTSmin defaults may-write to the combined matrix anyway (`pins.c`), which is a
  superset of the writes, so the pruning was never starved.
* The property is genuinely stutter invariant, so partial order reduction is
  legitimately applicable and the `-p` ITS-Tools passes is not the mistake. Spot
  on the source formula: `properties: stutter-invariant very-weak weak
  inherently-weak`.

The one asymmetry that is certain: LTSmin's own LTL front end refuses the
combination outright,

```
$B/pins2lts-seq-linux64 ./gal.so -p --pins-guards --when \
   --ltl='[]((<>(LTLAPp0==true)) || (X (X (<>(LTLAPp1==true)))))' --ltl-semantics=spin
** error **: The neXt operator is not allowed in combination with --por
```

while the same front end without `-p` finds the accepting cycle. That guard
(`pins2pins-ltl.c`) is syntactic -- it refuses any X, though this formula is
stutter invariant -- and it exists only on the path where LTSmin parses the
formula. Through `--hoa` LTSmin sees an automaton, the guard cannot fire, and
the reduction runs. `-m` also shows POR being set up before the automaton
exists: `Initializing POR dependencies: labels 109, guards 107`, our raw
counts, printed before `buchi has 3 states`.

Finding where the reduction actually loses the witness needs an instrumented
build -- print `ctx->visible` after `init_visible_labels` -- which
`~/git/LTSmin-BinaryBuilds` can produce: it clones upstream
`utwente-fmt/ltsmin` and already applies `patch/pins-impl.h`.

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
