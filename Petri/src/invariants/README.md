# `invariants/` — P/T flow and semi-flow solver

Computes a generative basis of the integer kernel of a sparse matrix (flows)
or of its non-negative part (semi-flows). Adapted from APT, heavily optimised.

* `InvariantMiddle.h` — front end: timeout wrapper, invariant printing,
  compression/decompression of permutation-equivalent invariants.
* `InvariantCalculator.h` — the elimination algorithm (Fourier–Motzkin style
  for semi-flows, Gaussian for flows) with pivot heuristics.
* `InvariantsTrivial.h` — culling of empty and duplicate columns before the
  main algorithm.
* `RowSigns.h`, `RowSignDomination.h` — per-row sign bookkeeping used to pick
  pivots and to detect dominated rows.
* `MixedSignsUniqueTable.h` — unique table for candidate rows.
* `Heuristic.h` — `EliminationHeuristic`: the option bundle for the solver.
