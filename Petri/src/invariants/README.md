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
* `Inequalities.h` — opt-in, read-only harvesting of one-sign certificates at
  phase-1 discard boundaries; see `algorithm.md`. No cone completeness claim.

`InvariantMiddle::computePInvariantsWithInequalities` returns `basis`,
`permutations`, and `inequalities` (`decreasing` and `increasing` matrices).
It takes the same variable-by-constraint matrix and heuristics as the regular
API, plus an optional cooperative deadline. For every returned nonnegative
column b, M^T b has the indicated sign. Coefficients use original variable
indices, independently of any basis compression. The legacy pair-returning
API remains unchanged and instantiates elimination without harvesting.

```cpp
auto result = petri::InvariantMiddle<long>::computePInvariantsWithInequalities(M, false);
// result.basis and result.permutations: ordinary flows.
// Each column b below uses M's original row indices, and satisfies M^T b <= 0.
for (const auto& b : result.inequalities.decreasing.getColumns()) {
    // Interpret b using the caller's variable names and initial marking.
}
```

The five nets in `../../examples/inequalities/` exercise named CLI reporting
with `--Pflows --collectInequalities` (or `--Psemiflows`).
