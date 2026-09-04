# `parse/` — input parsers

* `PTNetHandler.h` — expat SAX handler for PNML P/T nets (places, transitions,
  weighted arcs, initial markings; `toolspecific` is skipped as opaque).
* `PTNetLoader.h` — `loadXML(filename)`: drives expat over a file and returns
  a `SparsePetriNet<T>`.

* `NetResolver.h` — name-to-index lookup over a net and the desugaring of
  `is-fireable(t)` into `m[p] >= w` conditions, shared by the property
  parsers.
* `PropertyFile.h` — `loadPropertyFile(file, net, syntax)`: picks the parser
  by explicit syntax or by extension (`.xml` is MCC, anything else
  s-expressions).
* `mcc/` — MCC property XML parser producing `expr::Property` values
  (unsupported kinds kept with a comment).
* `sexpr/` — s-expression reader and property reader, the tool-to-tool
  syntax (`INTEROP.md` section 4).

A parser for another property syntax (the PetriVizu infix syntax) goes in its
own subfolder and must produce the same `expr::Property` values; the AST never
depends on the syntax it came from.
