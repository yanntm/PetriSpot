# Three total examinations

A note for ITS-Tools: add `UpperBoundsAll`, `StableMarkingAll` and
`QuasiLivenessAll`. Not implemented.

## Why

The MCC asks sixteen random conjunctions per model, so a score mixes model
hardness with formula luck. These three ask a question per place or per
transition, so the query set is determined by the model: what comes back is a
completion fraction that means something on its own, and comparable across
models and across the instances of a parametric family.

Three properties make them better instruments than the contest queries.

* **Half of each is self-certifying.** A quasi-live transition, an unstable
  place and a lower bound are all witnessed by a firing sequence; their
  negations need an argument. So the easy half needs no oracle at all, and the
  residue is by construction the part where search is not enough. The ratio
  says where PetriSpot can move the needle and where it is a spectator.
* **UpperBounds gives a distance.** `[lo, hi]` narrows even when nothing
  closes, so a heuristic can be tuned against a model it never solves. Boolean
  queries tell you nothing until they flip.
* **They land in the bulk regime.** Tens of thousands of goals on a large net
  is exactly the regime `PORTFOLIO.md` says changes how the portfolio must
  work; the experiment and the design stress the same thing.

## What they are

| examination | one query per | verdict | witnessed side |
| --- | --- | --- | --- |
| `QuasiLivenessAll` | transition *t* | can *t* fire in some reachable state | TRUE |
| `StableMarkingAll` | place *p* | does *p* hold the same marking in every reachable state | FALSE |
| `UpperBoundsAll` | place *p* | max of *p* over reachable markings | the lower end |

## How they run

The same path as the global properties, in `GlobalPropertySolver`: build the
cohort, discard what is cheap to discard, hand the rest to the solving loop.
The builders already exist in miniature — `buildQuasiLivenessProperty` emits one
`EF ENABLED(t)` per non-dominated transition, `buildStableMarkingProperty` one
per place — and the difference is only that the roll-up to a single examination
verdict is dropped: no early break, every atom is reported.

Cheap discards, before the loop:

* `QuasiLivenessAll` — `computeDominatedTransitions`; on a coloured model the
  skeleton settles negatives (dead in the skeleton is dead in the unfolding);
  structural and SMT dead-transition tests.
* `StableMarkingAll` — constant places are TRUE by construction; any place a
  first walk changes is FALSE by witness.
* `UpperBoundsAll` — constant places close at once; invariants give the upper
  side; NUPN or one-safe gives `≤ 1`.

## Output

One line per atom, MCC shape, named by the place or transition:

```
FORMULA <model>-QuasiLivenessAll-<transition> TRUE TECHNIQUES TOPOLOGICAL RANDOM_WALK
FORMULA <model>-StableMarkingAll-<place> FALSE TECHNIQUES ...
FORMULA <model>-UpperBoundsAll-<place> 12 TECHNIQUES ...
```

Verbose on a large net, and that is acceptable — the point of the exercise is
the whole vector. Two requirements on the output, without which the completion
fraction cannot be computed:

* **unanswered atoms are reported as such**, one line each, the way
  `--printUnknown` does in PetriSpot. Silence must not be readable as FALSE.
* **an unclosed bound reports its interval**, `lo` and `hi`, not nothing. That
  is the measurement the examination exists for.

Transition names are mangled on large nets, where composing them would blow up
in size; fall back to the index there and say so in the header line.

## What to measure

* Completion fraction per model, and the curve of atoms closed against time —
  a knee says a few pathological atoms, a slope says a wall.
* Witnessed against proved, per model.
* Whether the residue is the *same* atoms across the instances of a family. If
  three transitions resist in every `ResIsolation-N*P*`, that is structure, not
  a search failure.
* Where completion falls off a cliff along a parametric family.

## Caveat

There is no consensus oracle for a total query set. Correctness rests on the
self-certifying half, on the overlap with the sixteen contest formulas where
one exists, and on cross-tool agreement if another tool ever answers these. A
completion fraction is a measure of what we settled, never of what we settled
correctly, and should be reported as such.
