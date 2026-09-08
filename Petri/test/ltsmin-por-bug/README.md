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

## The suspect

Every reduced run prints

```
Visible groups: 0 / 111, labels: 2 / 111
```

Two state labels are visible -- the atomic propositions the automaton reads --
and no transition group is. A transition that changes an atomic proposition
must be visible or the reduction may drop the interleaving that exposes it, so
a visibility that never reaches the groups would let POR prune exactly the
transitions the property observes. That is a reading of the output, not a
diagnosis of LTSmin's code.

## What to do with it

Ours to work around, upstream to fix: the check is in
`LTSminRunner.checkProperty` (`if (doPOR && isStutterInvariant(pbody))` adds
`-p`), so dropping `-p`, or trusting a witness over an emptiness claim, are
both available. The reproduction folder built from the MCC archive is
`bench/models/Stigmergy-f03/` (git ignored), made by
`Petri/test/extract_property.py`; `Petri/test/logs/stig-dbg-1.log` is a full
run that takes the wrong branch, `Petri/test/logs/stigmergy-f03-debug.log` one
that answers FALSE by `STUTTER_TEST` before LTSmin is reached.
