# Standalone reduction command

`petri64 reduce -i model.pnml --props formulas.xml --output reduced`
reads the model and supported formulas, reduces their joint support, and writes
`reduced/model.pnet` and `reduced/properties.sexpr`. It performs no model checking.
`names.sexpr` retains the output-index to descriptive-name correspondence.

`--net` accepts an existing PNET instead of PNML. `--query N` selects one formula
before computing support; otherwise the whole file shares a reduction.
`--propsSyntax`, `--reductionMs`, and `--reductionNoAgglo` control input and effort.
The output directory must not already exist; its parent must exist. Write into
a staging directory and rename it after all files close successfully, so a
successful publication is a matched model/formula pair.

Solve later using any consumer of the existing exchange formats, for example:
`petri64 --net reduced/model.pnet --props reduced/properties.sexpr --totalTime=10`.
The output claims property preservation, not original state/arc counts or
original witness reconstruction; input counting records are discarded with a
diagnostic, since the rules these goals run cannot maintain them.

`petri64 reduce -i model.pnml --goal STATESPACE --output reduced` takes no
formulas: the net is reduced by the rules that keep the four StateSpace values
recoverable and `model.pnet` carries the counting record as named blocks
(`io/PNET.md`: `TMULT`, `PDROP`, `PCOEF`, `PCONST`), built by `CountingBlocks.h` from
`Counting.h`. A PNML input vouches for its own arcs (identity record); a
`--net` input is trusted only for the blocks it carries, so a PNET without
`TMULT` yields no `TMULT`. Unknown input blocks are dropped and named on
stderr. `hsc-pn --net reduced/model.pnet --states` then answers the
examination; `Petri/test/reduction/check_statespace.sh MODEL` runs that chain
against the deployed oracle.

`PropertyResults.h` is the analysis-side consumer of simplified properties:
scan for boolean bodies and exact bounds, print those results, and retain only
unresolved properties before constructing a solver. The standalone command
keeps simplified properties in its output instead of consuming them.
