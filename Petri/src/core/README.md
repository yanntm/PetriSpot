# `core/` — sparse substrate

Header-only templates on the integer type `T` (int, long, long long for the
three binaries). Everything else in the tree builds on these.

* `SparseArray.h` — sorted sparse vector (keys, values): merge-based `sumProd`,
  `greaterOrEqual` (galloping), `countContainsPos`, `manhattanDistance`,
  hashing and equality.
* `SparseBoolArray.h` — sorted sparse set of indices.
* `MatrixCol.h` — matrix stored by sparse columns; `transpose`, `sumProd`.
* `SparsePetriNet.h` — a P/T net: names, initial marking, `flowPT` and
  `flowTP` column matrices (one column per transition).
* `Arithmetic.hpp` — overflow-checked arithmetic and 128-bit printing.
* `Rational.h` — small rational type used by the invariant solver.
* `InvariantHelpers.h` — sparse vector kernels (`sumProdInto`, dot product).
* `Log.h` — timestamped INFO logging to standard output.
