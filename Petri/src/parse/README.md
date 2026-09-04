# `parse/` — input parsers

* `PTNetHandler.h` — expat SAX handler for PNML P/T nets (places, transitions,
  weighted arcs, initial markings; `toolspecific` is skipped as opaque).
* `PTNetLoader.h` — `loadXML(filename)`: drives expat over a file and returns
  a `SparsePetriNet<T>`.

* `mcc/` — MCC property XML parser producing `expr::Property` values
  (`is-fireable` desugared into input-arc conditions, unsupported kinds kept
  with a comment).

A parser for another property syntax (s-expressions, the PetriVizu infix
syntax) goes in its own subfolder and must produce the same `expr::Property`
values; the AST never depends on the syntax it came from.
