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
Input counting records cannot yet be maintained through this command; they are
explicitly discarded with a diagnostic. The output claims property preservation,
not original state/arc counts or original witness reconstruction.

`PropertyResults.h` is the analysis-side consumer of simplified properties:
scan for boolean bodies and exact bounds, print those results, and retain only
unresolved properties before constructing a solver. The standalone command
keeps simplified properties in its output instead of consuming them.
