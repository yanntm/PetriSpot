# `parse/` — input parsers

* `PTNetHandler.h` — expat SAX handler for PNML P/T nets (places, transitions,
  weighted arcs, initial markings; `toolspecific` is skipped as opaque).
* `PTNetLoader.h` — `loadXML(filename)`: drives expat over a file and returns
  a `SparsePetriNet<T>`.

Property parsers (MCC XML, s-expressions) live in subfolders as they arrive.
