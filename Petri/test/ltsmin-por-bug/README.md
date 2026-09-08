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

## Where the unsoundness is

Every reduced run prints

```
Visible groups: 0 / 111, labels: 2 / 111
```

LTSmin finds the two visible labels -- the atomic propositions the automaton
reads -- and no visible transition group. A transition that changes an atomic
proposition must be visible, or the reduction may drop the interleaving that
exposes it; with nothing visible the reduction is free to prune the witness.

The matrices we emit do carry the information. `label_matrix` says label 107
(`LTLAPp0`) reads state variables 12 and 28 and label 108 (`LTLAPp1`) reads 25;
`write_matrix` has 22 groups writing one of those three:

```
[1, 2, 3, 6, 7, 28, 29, 30, 31, 32, 40, 46, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63]
```

So the group visibility LTSmin reports as empty is derivable from what the model
declares. Two checks rule out the obvious explanations on our side:

* We register `GBsetDMInfoRead`, `GBsetDMInfoMustWrite` and `GBsetDMInfo`, never
  `GBsetDMInfoMayWrite` (`Gal2PinsTransformerNext`, `PetriNet2PinsTransformer`).
  Adding `GBsetDMInfoMayWrite(m, wm)` to `model.c` and rebuilding changes
  nothing: still `Visible groups: 0`, still `Empty product`.
* The property is genuinely stutter invariant, so partial order reduction is
  legitimately applicable and the `-p` ITS-Tools passes is not the mistake.
  Spot on the source formula: `properties: stutter-invariant very-weak weak
  inherently-weak`.

LTSmin's own LTL front end refuses the combination outright:

```
$B/pins2lts-seq-linux64 ./gal.so -p --pins-guards --when \
   --ltl='[]((<>(LTLAPp0==true)) || (X (X (<>(LTLAPp1==true)))))' --ltl-semantics=spin
** error **: The neXt operator is not allowed in combination with --por
```

and without `-p` that same front end finds the accepting cycle. The guard is
syntactic -- it refuses any X, though this formula is stutter invariant -- and
it only exists on the path where LTSmin parses the formula. Through `--hoa` it
sees an automaton, the guard cannot fire, and the reduction runs with no visible
group.

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
