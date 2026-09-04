# `io/` — exporters

* `SparseMatrixIO.h` — KERS binary sparse matrix format (see `KERS.md`); file
  and stream overloads, the stream ones read or write one block.
* `PNETIO.h` — PNET binary net: header plus three KERS blocks (flowPT,
  flowTP, marking); the tool-to-tool net format of `INTEROP.md`.
* `MatrixExporter.h` — ASCII sparse matrix export.
* `PNMLExport.h` — normalised PNML export of a net.
* `FlowPrinter.h` — Graphviz dot rendering of a net with highlighted flows.
