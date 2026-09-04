/*
 * InvariantDriver.h
 *
 * Invariant computations from the command line: over a KERS matrix (no net)
 * or over a loaded net, with listing and KERS export of the basis.
 */
#ifndef PETRI_CLI_INVARIANTDRIVER_H_
#define PETRI_CLI_INVARIANTDRIVER_H_

#include <chrono>
#include <iostream>
#include <string>
#include <vector>

#include "cli/Options.h"
#include "core/MatrixCol.h"
#include "core/SparsePetriNet.h"
#include "invariants/InvariantMiddle.h"
#include "io/SparseMatrixIO.h"

namespace petri::cli
{

inline long millisSince (std::chrono::steady_clock::time_point t)
{
  return std::chrono::duration_cast<std::chrono::milliseconds> (std::chrono::steady_clock::now () - t).count ();
}

/**
 * --loadKERS path: the matrix is the incidence matrix, rows are variables.
 * Returns the process exit code.
 */
template<typename T>
  int runInvariantsFromKERS (const Options &o)
  {
    if (!o.invariants ()) {
      std::cerr << "Error: --loadKERS requires a flow/semiflow flag (--Pflows, --Psemiflows, --Tflows, --Tsemiflows)."
          << std::endl;
      return 1;
    }
    if (o.useCompression) {
      std::cerr << "Warning: --useCompression is not supported with --loadKERS (no net structure); flag ignored."
          << std::endl;
    }
    MatrixCol<T> M = SparseMatrixIO<T>::read (o.loadKERSFile);
    if (M.getColumnCount () == 0 && M.getRowCount () == 0) return 1;
    bool transposed = o.tflows || o.tsemiflows;
    bool semi = o.psemiflows || o.tsemiflows;
    MatrixCol<T> toUse = transposed ? M.transpose () : M;
    EliminationHeuristic heur = o.heuristic (false);
    auto time = std::chrono::steady_clock::now ();
    auto [mat, perms] = InvariantMiddle<T>::computePInvariants (toUse, semi, o.timeout, heur);
    std::cout << "Computed " << mat.getColumnCount () << " " << (transposed ? "T" : "P") << " "
        << (semi ? "semi" : "") << "flows in " << millisSince (time) << " ms." << std::endl;
    // basisKERS implies program-to-program mode: no human-readable listing
    if (!o.quiet && o.basisKERSFile.empty ()) {
      size_t ndim = M.getRowCount ();
      std::vector<std::string> names;
      names.reserve (ndim);
      for (size_t i = 0; i < ndim; ++i) names.push_back ("x" + std::to_string (i));
      std::vector<T> zeros (ndim, 0);
      InvariantMiddle<T>::printInvariant (mat, perms, names, zeros);
    }
    if (!o.basisKERSFile.empty ()) {
      if (!SparseMatrixIO<T>::write (mat, o.basisKERSFile)) return 1;
      std::cout << "Exported basis to " << o.basisKERSFile << std::endl;
    }
    return 0;
  }

/** One flow kind on a net: compute, report, list or export. */
template<typename T>
  void computeFlows (const Options &o, const SparsePetriNet<T> &pn, bool onTransitions, bool semi)
  {
    EliminationHeuristic heur = o.heuristic (o.useCompression);
    auto time = std::chrono::steady_clock::now ();
    MatrixCol<T> sumMatrix = MatrixCol<T>::sumProd (-1, pn.getFlowPT (), 1, pn.getFlowTP ());
    if (onTransitions) sumMatrix = sumMatrix.transpose ();
    auto [mat, perms] = InvariantMiddle<T>::computePInvariants (sumMatrix, semi, o.timeout, heur);
    std::cout << "Computed " << mat.getColumnCount () << " " << (onTransitions ? "T " : "P ")
        << (semi ? "semi" : "") << "flows ";
    if (!perms.empty ()) {
      std::cout << "with " << perms.size () << " permutations.\n";
      std::cout << "Total decompressed invariants :" << InvariantMiddle<T>::countInvariant (mat, perms) << std::endl;
    }
    std::cout << " in " << millisSince (time) << " ms." << std::endl;
    if (!o.quiet && o.basisKERSFile.empty ()) {
      if (onTransitions) {
        std::vector<T> zeros (pn.getTnames ().size (), 0);
        InvariantMiddle<T>::printInvariant (mat, perms, pn.getTnames (), zeros);
      } else {
        InvariantMiddle<T>::printInvariant (mat, perms, pn.getPnames (), pn.getMarks ());
      }
    }
    if (!o.basisKERSFile.empty ()) {
      SparseMatrixIO<T>::write (mat, o.basisKERSFile);
      std::cout << "Exported basis to " << o.basisKERSFile << std::endl;
    }
  }

/** Every flow kind requested by the options, P side then T side. */
template<typename T>
  void runInvariants (const Options &o, const SparsePetriNet<T> &pn)
  {
    if (o.pflows || o.psemiflows) computeFlows (o, pn, false, o.psemiflows);
    if (o.tflows || o.tsemiflows) computeFlows (o, pn, true, o.tsemiflows);
  }

} // namespace petri::cli

#endif /* PETRI_CLI_INVARIANTDRIVER_H_ */
