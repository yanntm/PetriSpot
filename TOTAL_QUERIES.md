# Three total examinations

ITS-Tools examinations `QuasiLivenessAll`, `StableMarkingAll` and
`UpperBoundsAll`, P/T nets only. Implemented in
`fr.lip6.move.gal.application.solver.total` (`TotalExamination`,
`TotalSolver`, `TotalPrinter`); the atoms come from
`solver.global.GlobalAtoms`, shared with the roll-up examinations.

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

As an MCC examination: `-examination QuasiLivenessAll` builds one atom per
object on the net as parsed, in place of reading a property file, and falls
through to the `ReachabilityFireability`, `ReachabilityCardinality` or
`UpperBounds` pipeline with all its engines. Nothing is discarded up front:
the dominated-transition mask of the roll-up path is an implication, not an
answer, and the pipeline settles constant places and dead transitions on its
own. The `-timeout` flag is a per-engine budget, not a deadline; the run
ends when the cohort is closed or when the harness kills it.

## Output

One header, then one line per atom the moment its verdict lands. Objects
are named by definition index in the PNML, `p<i>` and `t<i>`; names can be
resolved from the file afterwards if ever needed.

```
TOTAL QuasiLivenessAll 1242
QLIVE t17 TRUE TOPOLOGICAL RANDOM_WALK
STABLE p3 FALSE DECISION_DIAGRAMS TOPOLOGICAL
BOUND p9 12 TOPOLOGICAL SAT_SMT
BOUND p11 ? 5 inf
```

The last shape is an open bound, printed at every checkpoint of the bounds
solver where its `[max seen, structural bound]` pair moved; a reader keeps
the last line of each object. Nothing is printed post hoc, so an external
kill loses at most the current phase. An atom with no line is unanswered:
the header gives the count the completion fraction is taken against.

For the cluster, `MCC-drivers/gen_total_oracles.sh` writes one oracle per
examination and instance with a `?` per object (`<instance>-QLA.out`,
`-SMA.out`, `-UBA.out`), so `run_test.pl` runs them like any examination and
reports the unanswered atoms as missing results. A `?` is to be replaced by
a verdict once one is trusted.

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
