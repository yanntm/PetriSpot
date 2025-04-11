#ifndef INVARIANTCALCULATOR_H_
#define INVARIANTCALCULATOR_H_

/**
 * A calculator for invariants computing flows or semiflows.
 * Provides a sparse implementation of a variant on Farkas algorithm.
 *
 * The algorithm proceed in two phases:
 * - Phase 1 computes a basis of the kernel of the input matrix, which is a basis of the set of flows.
 * This is worst case proportional to the number of columns of the matrix.
 * - Phase 2 massages the flows into semiflows. This is worst case exponential in the number of flows.
 *
 * Originally, this file is derived from an implementation provided by the APT tool,
 * from the team of Eike Best, University of Oldenburg, in Java.
 * Dennis-Michael Borde, Manuel Gieseking are credited as authors in the source code.
 *
 * The algorithm itself is based on the PIPE algorithm described in
 * http://pipe2.sourceforge.net/documents/PIPE-Report.pdf (page 19) which is
 * based on the paper of D'Anna and Trigila "Concurrent system analysis using
 * Petri nets – an optimised algorithm for finding net invariants", Mario D'Anna
 * and Sebastiano Trigila, Computer Communications vol 11, no. 4 august 1988.
 * according to source code comments.
 *
 * It was first heavily adapted for ITS-Tools by Yann Thierry-Mieg starting from 2017,
 * leading to a sparse implementation of the algorithm described in the PNSE23 paper
 * "Efficient Strategies to Compute Invariants, Bounds and Stable Places of Petri nets"
 * by Yann Thierry-Mieg.
 *
 * In 2024, the code was ported to C++ (with the help of Soufiane El Mahdi M1 student at Sorbonne University).
 *
 * In beginning of 2025 the code was practically entirely rewritten, though main elements of the algorithm remain, the source
 * is barely recognizable. At this occasion several new heuristics and performance enhancements were added, implementing
 * new ideas developed with the help of Denis Poitrenaud and Emmanuel Paviot-Adet.
 *
 * @author Yann Thierry-Mieg.
 *
 */

#include <numeric>
#include <stdexcept> // for std::overflow_error
#include <string>
#include <iostream>
#include <vector>
#include <unordered_set>
#include "SparseArray.h"
#include "MatrixCol.h"
#include "SparseBoolArray.h"
#include "Arithmetic.hpp"
#include "InvariantHelpers.h"
#include "RowSigns.h"
#include "RowSignDomination.h"
#include "InvariantsTrivial.h"
#include "Heuristic.h"
#include "MixedSignsUniqueTable.h"

namespace petri
{

/**
 * @class InvariantCalculator
 * @brief Computes invariants (flows or semiflows) of a Petri net using a sparse matrix representation.
 *
 * This class provides a static interface to calculate invariants from a given incidence matrix. It supports both flows
 * (allowing negative coefficients) and semiflows (non-negative coefficients only). The algorithm operates in two phases:
 * - Phase 1: Reduces the matrix to a basis of invariants using row elimination.
 * - Phase 2: Optionally refines the basis into a minimal set of semiflows with positive coefficients.
 *
 * The implementation uses sparse data structures (SparseArray, MatrixCol) and incorporates heuristics for pivot selection
 * and row elimination to optimize performance.
 *
 * @tparam T The numeric type for matrix entries (e.g., int, long).
 */
template<typename T>
  class InvariantCalculator
  {
    // Enable for verbose debugging output during development.
    static inline const bool DEBUG = false;

    /**
     * Hidden constructor
     */
    InvariantCalculator ()
    {
    }
  public:

    // Permutation struct for symmetry axes
    struct Permutation
    {
      size_t index;                    // New row index (e.g., perm_42)
      std::vector<SparseArray<T>> elements; // Positive terms (e.g., {p1}, {2p3 + p4})

      Permutation (size_t idx, const std::vector<SparseArray<T>> &pv)
          : index (idx), elements (pv)
      {
      }

      friend std::ostream& operator<< (std::ostream &os, const Permutation &p)
      {
        os << "Permutation " << p.index << " : ";
        for (const auto &el : p.elements) {
          os << el << " ";
        }
        return os;
      }
    };
    using Permutations = std::vector<Permutation>;

    /**
     * @brief Computes a generating set of invariants (flows or semiflows) from a Petri net incidence matrix.
     *
     * This method implements a sparse variant of the PIPE algorithm.
     * It operates in two phases:
     * 1. Phase 1: Reduces the input matrix to a basis of flows using row elimination (see phase1PIPE).
     * 2. Phase 2: If onlyPositive is true, refines the basis into a minimal set of semiflows (see phase2Pipe).
     *
     * The input matrix is normalized (removing duplicate and scaled columns) before processing. The result is a set of sparse
     * vectors representing the invariants.
     *
     * @param mat          The incidence matrix (rows = places, columns = transitions).
     * @param onlyPositive If true, computes semiflows (non-negative coefficients); if false, computes flows (allowing negatives).
     * @param heur         Heuristic settings for row elimination and pivot selection (default: EliminationHeuristic()).
     * @return An unordered set of SparseArray<T> representing the invariants.
     */
    static std::pair<MatrixCol<T>, Permutations> calcInvariantsPIPE (
        MatrixCol<T> &mat, bool onlyPositive, const EliminationHeuristic &heur =
            EliminationHeuristic ())
    {
      if (mat.getColumnCount () == 0 || mat.getRowCount () == 0) {
        return {MatrixCol<T> (), {}};
      }

      // while it's not classical, we actually work with columns in the algorithm
      // The matrix data structure is asymetric in complexity between rows and columns
      mat = mat.transpose ();

      size_t colCount = mat.getColumnCount ();
      // Normalize the transposed matrix: remove duplicate columns and scale them to canonical form.
      mat.normalizeAndReduce (false);

      if (mat.getColumnCount () < colCount) {
        std::cout << "Normalized transition count is " << mat.getColumnCount ()
            << " out of " << colCount << " initially." << std::endl;
      }

      mat = mat.transpose ();

      MatrixCol<T> matB = phase1PIPE (mat, onlyPositive, heur);

//		const MatrixCol<T> matB = phase1PIPE(new MatrixCol<T>(mat));
      // We want to work with columns in this part of the algorithm
      // We add and remove columns all day => we want to switch to a column based
      // representation
      // order of rows is really irrelevant + columns which are identical up to
      // scaling factor are useless
      // let's use a set of columns.
      matB.normalizeAndReduce (true);

      if (!onlyPositive) {
        return {matB, {}};
      }

      // phase 2
      std::cout << "// Phase 2 : computing semi flows from basis of "
          << matB.getColumnCount () << " invariants " << std::endl;

      /* FACTORIZATION */
      if (heur.useCompression ()) {
        /* FACTORIZATION */
        auto [perms, colsB] = factorizeBasis (matB);
        matB = colsB;
        phase2Pipe (matB, heur);
        return {matB,perms};
      } else {
        phase2Pipe (matB, heur);
        return {matB, {}};
      }

    }

  private:

    /**
     * @brief Phase 1 of the invariant calculation: computes a basis of flows by eliminating rows.
     *
     * This method reduces the input matrix `matC` to a zero matrix while updating a transformation matrix `matB` to
     * represent a basis of flows. It uses row elimination with pivot selection guided by heuristics. Key steps:
     * - Removes trivial (empty) and duplicate columns early.
     * - Discards rows that cannot contribute to semiflows (if onlyPositive is true).
     * - Iteratively eliminates rows until matC is zero, tracking changes in matB.
     *
     * @param matC         The incidence matrix to reduce (modified in place).
     * @param onlyPositive If true, enforces constraints for semiflows during elimination.
     * @param heur         Heuristic settings for pivot selection and elimination order.
     * @return The transformation matrix matB representing the basis of flows.
     */
    static MatrixCol<T> phase1PIPE (MatrixCol<T> &matC, bool onlyPositive,
                                    const EliminationHeuristic &heur)
    {
      // Build the initial transformation matrix.
      MatrixCol<T> matB = MatrixCol<T>::identity (matC.getColumnCount (),
                                                  matC.getColumnCount ());

      // Vector to hold trivial invariants.
      std::vector<SparseArray<T>> trivialInv;

      if (heur.useCulling ()) {
        auto time = std::chrono::steady_clock::now ();

        // Remove trivial invariants (empty columns in matC) early.
        size_t culled = cullConstantColumns (matC, matB, trivialInv);
        // Remove duplicate columns
        culled += cullDuplicateColumns (matC, matB, trivialInv);

        std::cout << "Culled " << culled << " columns in "
            << std::chrono::duration_cast < std::chrono::milliseconds
            > (std::chrono::steady_clock::now () - time).count () << " ms."
                << std::endl;
      } else {
        std::cout << "Skipping cull step" << std::endl;
      }

      std::cout << "// Phase 1: matrix " << matC.getRowCount () << " rows "
          << matC.getColumnCount () << " cols " << matC.getEntryCount ()
          << " entries" << std::endl;
      RowSigns rowSigns (matC, heur.useSingleSignRow ());

      // For semiflows, discard rows with all-positive or all-negative entries early,
      // as they cannot balance to zero with non-negative coefficients.
      if (onlyPositive) {
        size_t elim = 0;
        for (size_t row = 0; row < matC.getRowCount (); row++) {
          const auto &rs = rowSigns.get (row);
          auto pm = rs.pMinus.size ();
          auto pp = rs.pPlus.size ();
          if ((pm == 0 && pp > 0) || (pp == 0 && pm > 0)) {
            eliminateRowWithPivot (PivotChoice (row, -1), matC, matB, rowSigns,
                                   onlyPositive);
            elim++;
          }
        }
        if (elim) {
          std::cout << "// After discarding improbable semiflows : matrix "
              << matC.getRowCount () << " rows " << matC.getColumnCount ()
              << " cols " << matC.getEntryCount () << " entries" << std::endl;
        }
      }

      // stats on which rule applied
      std::pair<size_t, size_t> counts (0, 0);

      // Main elimination loop: reduce matC to zero by applying row operations, updating matB accordingly.
      while (!matC.isZero ()) {
        applyRowElimination (matC, matB, rowSigns, counts, onlyPositive, heur);
        if (DEBUG) {
          std::cout << "Mat max : " << matC.maxVal () << std::endl;
          std::cout << "B max : " << matB.maxVal () << std::endl;
        }
      }
      std::cout << "Finished phase 1 with " << counts.first
          << " SingleSign rule and " << counts.second << " generalized "
          << std::endl;
      // Re-add the trivial invariants back into matB.
      for (const auto &inv : trivialInv) {
        matB.appendColumn (inv);
      }

      return matB;
    }

    /**
     * @brief Phase 2 of the invariant calculation: refines a basis of flows into semiflows.
     *
     * This method takes a basis of flows (colsB) and computes a minimal set of semiflows (non-negative coefficients) by:
     * - Clearing rows with only negative entries (and their associated columns).
     * - Iteratively eliminating rows using linear combinations to ensure all coefficients are positive.
     * - Normalizing the resulting vectors.
     *
     * @param colsB       The basis of flows from Phase 1 (modified in place).
     * @param semiFlows   Output set of semiflows.
     * @param heur        Heuristic settings for row selection.
     */
    static void phase2Pipe (MatrixCol<T> &colsB,
                            const EliminationHeuristic &heur)
    {
      RowSignsDomination<T> rowSigns (colsB);
      std::unordered_set<size_t> basisIndices;
      MixedSignsUniqueTable<T> msut (colsB); // Initialize MSUT regardless, but only seed in Q+
      size_t filteredCount = 0;

      // Step 1: Remove rows with only negative entries, as they cannot contribute to semiflows.
      {
        std::unordered_set<size_t> tokill;
        for (const auto &rs : rowSigns) {
          if (rs.pPlus.size () == 0) {
            for (size_t i = 0, ie = rs.pMinus.size (); i < ie; i++) {
              tokill.insert (rs.pMinus.keyAt (i));
            }
          }
        }
        if (tokill.size () > 0) {
          for (const auto &col : tokill) {
            clearColumn (col, colsB, rowSigns);
          }
          if (DEBUG) {
            std::cout << "Cleared " << tokill.size ()
                << " cols due to row with only negative entries\n";
          }
        }
      }

      // Step 2: Seed basisIndices with all-positive vectors and msut with mixed-sign vectors (Q+ only)
      if (heur.useQPlusBasis ()) {
        std::cout << "Seeding basis and mixed-sign table with "
            << colsB.getColumnCount () << " columns\n";
        for (size_t i = 0; i < colsB.getColumnCount (); ++i) {
          if (colsB.getColumn (i).isPurePositive ()) {
            basisIndices.insert (i);
          }
          msut.insert (i);
        }
        if (DEBUG) {
          std::cout << "Initial basis size: " << basisIndices.size () << "\n";
          for (const auto &index : basisIndices) {
            std::cout << index << " :" << colsB.getColumn (index) << "\n";
          }
          std::cout << "Initial mixed-sign unique table size: " << msut.size ()
              << "\n";
          std::cout << "End of initial seeding" << std::endl;
        }
      }

      std::unordered_set<size_t> posRows;
      for (size_t row = 0; row < colsB.getRowCount (); row++) {
        if (rowSigns.get (row).pMinus.size () == 0) {
          posRows.insert (row);
        }
      }

      int iter = 0;
      // Step 3: Iteratively eliminate rows to ensure all coefficients are non-negative.
      while (true) {
        ssize_t tRow = -1;

        if (heur.useSingleSignRow ()) {
          rowSigns.findSingleSignRow (heur.getLoopLimit ());
        }

        if (tRow == -1) {
          // look for a "small" row
          tRow = findBestFMERow (colsB, rowSigns, heur.getLoopLimit ());
        }
        if (tRow == -1) {
          break;
        }
        eliminateRowFME (tRow, colsB, rowSigns, basisIndices, heur,
                         filteredCount, msut, posRows);  // Pass msut always

        if (DEBUG && ++iter /*(++iter % 10 == 0)*/) {
          std::cout << "Iteration " << iter << "Basis size: "
              << basisIndices.size () << "\n";
          for (const auto &index : basisIndices) {
            std::cout << index
                << " " /*" :" << colsB.getColumn (index) << "\n"*/;
          }
          std::cout << "\nMixed-sign unique table size: " << msut.size ()
              << "\n";
          std::cout << "Mixed-sign stats - attempted insertions: "
              << msut.getAttemptedInsertions () << ", successful insertions: "
              << msut.getSuccessfulInsertions () << "\n";
          RowSignsDomination<T>::confrontIndex (colsB, rowSigns, basisIndices);
        }
      }

      // Step 4: Stats and finalize
      std::cout << "After removing all negative signs "
          << colsB.getColumnCount () << " entries in matrix.\n";
      if (heur.useQPlusBasis ()) {
        std::cout << "Q+ minimal vectors: " << basisIndices.size ()
            << ", filtered " << filteredCount << " non-minimal vectors\n";
        std::cout << "Mixed-sign stats - attempted insertions: "
            << msut.getAttemptedInsertions () << ", successful insertions: "
            << msut.getSuccessfulInsertions () << "\n";

        colsB.dropEmptyColumns ();
        std::cout << "After dropping empty entries : "
            << colsB.getColumnCount () << " semiflows\n";
        if (DEBUG) {
	        colsB = minimizeBasisWithSupport(colsB);
        	std::cout << "After minimization with support: "
                    << colsB.getColumnCount () << " semiflows\n";
	}
      } else if (heur.useMinimization ()) {
        colsB = minimizeBasisWithSupport(colsB);
      } else {
        colsB.normalizeAndReduce (true);
      }
      std::cout << "Final semi flow basis size: " << colsB.getColumnCount ()
          << "\n";
    }

    static MatrixCol<T> minimizeBasisWithSupport (MatrixCol<T> &colsB)
    {
      auto startTime = std::chrono::high_resolution_clock::now ();

      // Step 0: Normalize input
      std::cout
          << "Normalizing (culling duplicates and empty) from input matrix with "
          << colsB.getColumnCount () << " columns";
      colsB.normalizeAndReduce ();
      std::cout << "produced a matrix with " << colsB.getColumnCount ()
          << " columns \n";

      // Step 1: Sort colsB by support size (descending)
      colsB.sortByColumnSize (true); // Largest to smallest, smallest at high indices

      // Step 2: Build row index
      RowSigns<T> rowIndex (colsB, false);  // All non-zero rows

      // Step 3: Process vectors
      MatrixCol<T> basis (colsB.getRowCount ());
      std::cout << "Starting minimization loop with " << colsB.getColumnCount ()
          << " columns\n";

      for (ssize_t lastVictim = colsB.getColumnCount () - 1; lastVictim >= 0;
          --lastVictim) {
        // Skip if column is empty
        if (colsB.getColumn (lastVictim).size () == 0) {
          continue;
        }

        if (DEBUG) {
          // Adopt vector
          std::cout << "Adopting column " << lastVictim << " with support size "
              << colsB.getColumn (lastVictim).size () << ": "
              << colsB.getColumn (lastVictim) << "\n";
        }
        basis.appendColumn (colsB.getColumn (lastVictim));  // Copy

        // Compute intersection
        const SparseArray<T> &currentCol = colsB.getColumn (lastVictim);
        SparseBoolArray intersectingCols =
            rowIndex.get (currentCol.keyAt (0)).pPlus;
        intersectingCols.put (lastVictim, false); // Remove self
        for (size_t i = 1; i < currentCol.size (); ++i) {
          size_t row = currentCol.keyAt (i);
          intersectingCols.restrict (rowIndex.get (row).pPlus);
        }
        if (DEBUG) {
          std::cout << "Intersection for column " << lastVictim << " contains "
              << intersectingCols.size () << " columns: ";
          for (size_t j = 0; j < intersectingCols.size (); ++j) {
            std::cout << intersectingCols.keyAt (j) << " ";
          }
          std::cout << "\n";
        }

        intersectingCols.put (lastVictim, true); // Readd self
        // Discard redundant vectors (including self) in decreasing order
        for (ssize_t j = intersectingCols.size () - 1; j >= 0; --j) {
          size_t colIdx = intersectingCols.keyAt (j);
          const SparseArray<T> &candidateCol = colsB.getColumn (colIdx);
          if (DEBUG && candidateCol.size () >= 0 && colIdx != lastVictim) {
            std::cout << "Eliminating column " << colIdx << " (" << candidateCol
                << ") because it is dominated by adopted column " << lastVictim
                << " (" << currentCol << ")\n";
          }
          InvariantCalculator<T>::clearColumn (colIdx, colsB, rowIndex);
        }
      }

      // Step 4: Finalize and return
      auto endTime = std::chrono::high_resolution_clock::now ();
      auto duration = std::chrono::duration_cast < std::chrono::milliseconds
          > (endTime - startTime);
      std::cout << "Minimization complete. Basis size: "
          << basis.getColumnCount () << ". Time: " << duration.count ()
          << " ms\n";

      return basis;
    }

    /**
     * @brief Minimizes a basis of positive semiflows into a minimal generating set with tracing.
     *
     * Processes colsB into a minimal basis in semiFlows, logging initial non-empty columns,
     * final insertions, duplicates removed, and execution time.
     *
     * @param colsB     Input basis of semiflows (positive, normalized).
     * @param semiFlows Output set of minimal semiflows.
     */
    // Hash functor for std::pair<size_t, size_t> to use in unordered_map
    struct PairHash
    {
      std::size_t operator() (const std::pair<size_t, size_t> &p) const
      {
        // Simple hash combination: shift first and XOR with second
        return (p.first << 16) ^ p.second;
      }
    };

    static MatrixCol<T> minimizeBasis (MatrixCol<T> &colsB)
    {
      auto startTime = std::chrono::steady_clock::now (); // Start timing

      // start with a normalization
      colsB.normalizeAndReduce (true);

      // Step 1: Collect non-empty columns into partitioned map
      std::unordered_map<std::pair<size_t, size_t>,
          std::vector<SparseArray<T>*>, PairHash> candidatesMap;
      size_t initialNonEmpty = 0; // Count non-empty columns (XX)
      for (size_t i = 0; i < colsB.getColumnCount (); i++) {
        SparseArray<T> *col = &colsB.getColumn (i);
        if (col->size () > 0) { // Skip empty
          std::pair<size_t, size_t> key =
            { col->keyAt (0), col->keyAt (col->size () - 1) };
          candidatesMap[key].push_back (col);
          initialNonEmpty++;
        }
      }

      // Step 2: Initialize visitation with all keys
      std::unordered_set<std::pair<size_t, size_t>, PairHash> indexesToBeVisited;
      for (const auto &pair : candidatesMap) {
        indexesToBeVisited.insert (pair.first);
      }

      // Compute runtime
      auto dur1 = std::chrono::duration_cast < std::chrono::milliseconds
          > (std::chrono::steady_clock::now () - startTime).count ();
      std::cout << "Preparing partition took :" << dur1 << " ms" << std::endl;
      if (DEBUG) {
        for (const auto &ent : candidatesMap) {
          std::cout << ent.second.size () << ", ";
        }
        std::cout << std::endl;
      }

      size_t iterCount = 0;
      // Step 3: Domination testing with fixpoint loop
      while (!indexesToBeVisited.empty ()) {
        iterCount++;
        // Sort keys by partition size for this pass
        std::vector<std::pair<size_t, size_t>> keys (
            indexesToBeVisited.begin (), indexesToBeVisited.end ());
        std::sort (
            keys.begin (), keys.end (),
            [&] (std::pair<size_t, size_t> a, std::pair<size_t, size_t> b) {
              return candidatesMap[a].size() < candidatesMap[b].size();
            });

        // Process each key
        for (std::pair<size_t, size_t> key : keys) {
          indexesToBeVisited.erase (key);
          auto &vecs = candidatesMap[key];
          std::vector<char> skips (vecs.size (), 0);

          std::sort (vecs.begin (), vecs.end (),
                     [] (const SparseArray<T> *a, const SparseArray<T> *b) {
                       return a->size() < b->size();
                     });

          for (size_t i = 0; i < vecs.size (); i++) {
            if (skips[i]) continue;
            for (size_t j = i + 1; j < vecs.size (); j++) {
              if (skips[j]) continue;
              if (subtractAllFrom (vecs[i], vecs[j])) {
                skips[j] = 1;
                continue;
              }
              if (subtractAllFrom (vecs[j], vecs[i])) {
                skips[i] = 1;
                break;
              }
            }
          }

          bool revisitCurrent = false;
          for (size_t k = vecs.size (); k-- > 0;) {
            if (!skips[k]) continue;
            SparseArray<T> *vec = vecs[k];
            if (vec->size () == 0) {
              // Discard empty
            } else if (vec->keyAt (0) != key.first
                || vec->keyAt (vec->size () - 1) != key.second) {
              std::pair<size_t, size_t> newKey =
                { vec->keyAt (0), vec->keyAt (vec->size () - 1) };
              candidatesMap[newKey].push_back (vec);
              indexesToBeVisited.insert (newKey);
            } else {
              revisitCurrent = true;
            }
            vecs.erase (vecs.begin () + k);
          }
          if (revisitCurrent && !vecs.empty ()) {
            indexesToBeVisited.insert (key);
          }
          if (vecs.empty ()) {
            candidatesMap.erase (key);
          }
        }
      }

      // Step 4: Transfer to semiFlows with tracing
      size_t insertions = 0; // Total insertions attempted (YY)
      MatrixCol<T> semiFlows (colsB.getRowCount ());
      for (const auto &pair : candidatesMap) {
        for (SparseArray<T> *vec : pair.second) {
          if (vec->size () > 0) {
            semiFlows.appendColumn (std::move (*vec)); // Attempt insertion
            insertions++;
          }
        }
      }

      size_t finalSize = semiFlows.getColumnCount (); // Final unique flows (ZZ)
      size_t duplicatesRemoved = insertions - finalSize; // Duplicates (YY - ZZ if no prior entries)

      // Compute runtime
      auto duration = std::chrono::duration_cast < std::chrono::milliseconds
          > (std::chrono::steady_clock::now () - startTime).count ();

      // Output trace message
      std::cout << "Minimization of " << initialNonEmpty
          << " non-empty flows to a minimal basis of " << finalSize << " flows"
          << (duplicatesRemoved > 0 ?
              " (" + std::to_string (duplicatesRemoved)
                  + " duplicates removed)" :
              "") << " took " << iterCount << " iterations " << " in "
          << duration << " ms" << std::endl;

      return semiFlows;
    }
    // Helper: Subtracts coli from colj as many times as possible
    static bool subtractAllFrom (const SparseArray<T> *coli,
                                 SparseArray<T> *colj)
    {
      T count = SparseArray<T>::countContainsPos (*colj, *coli);
      if (count == 0) return false;
      sumProdIntoNoChange<T> (1, *colj, -count, *coli); // colj -= count*coli
      normalize (*colj);
      return true;
    }

    /**
     * @brief Eliminates a target row in Phase 2 using Fourier Motzkin Elimination (FME) variant.
     *
     * Combines columns from P+ (positive entries) and P- (negative entries) in the target row to produce new columns
     * with non-negative coefficients. If a pure positive column exists, it is used as a pivot to avoid introducing
     * new negatives. The resulting columns are normalized and appended or updated in place.
     *
     * @param targetRow The row to eliminate.
     * @param colsB     The basis matrix (modified in place).
     * @param rowSigns  Bookkeeping of row signs (updated during elimination).
     */
    static void eliminateRowFME (ssize_t targetRow, MatrixCol<T> &colsB,
                                 RowSignsDomination<T> &rowSigns,
                                 std::unordered_set<size_t> &basisIndices,
                                 const EliminationHeuristic &heur,
                                 size_t &filteredCount,
                                 MixedSignsUniqueTable<T> &msut,
                                 std::unordered_set<size_t> &posRows)
    {
      const auto &rs = rowSigns.get (targetRow);
      if (DEBUG) {
        std::cout << "Current basis indexes " << basisIndices.size () << "\n";

        std::cout << "Eliminating row " << targetRow << " with "
            << rs.pPlus.size () << " plus and " << rs.pMinus.size ()
            << " minus columns\n";
        std::cout << rs << "\n";
        for (size_t i = 0, ie = rs.pPlus.size (); i < ie; i++) {
          std::cout << "P+ " << rs.pPlus.keyAt (i) << " : "
              << colsB.getColumn (rs.pPlus.keyAt (i)) << "\n";
        }
        for (size_t i = 0, ie = rs.pMinus.size (); i < ie; i++) {
          std::cout << "P- " << rs.pMinus.keyAt (i) << " : "
              << colsB.getColumn (rs.pMinus.keyAt (i)) << "\n";
          if (colsB.getColumn (rs.pMinus.keyAt (i)).size () == 0) {
            std::cout << "Empty column detected\n";
            exit (2);
          }
        }

        std::cout << std::flush;
      }

      if (rs.pPlus.size () == 0) {
        SparseBoolArray toVisit = rs.pMinus;
        for (size_t i = 0; i < toVisit.size (); ++i) {
          msut.erase (i);
          clearColumnWithBasis (i, colsB, rowSigns, false); // Non-basis by default
        }
        if (DEBUG) {
          std::cout << "Cleared row " << targetRow << "\n";
        }
        return;
      }

      ssize_t purePos = -1;
      if (!heur.useQPlusBasis () && rs.pPlus.size () > 1) {
        for (size_t i = 0; i < rs.pPlus.size (); ++i) {
          size_t candCol = rs.pPlus.keyAt (i);
          const auto &col = colsB.getColumn (candCol);
          if (col.isPurePositive ()) {
            if (purePos == -1
                || colsB.getColumn (purePos).size () > col.size ()) {
              purePos = candCol;
            }
          }
        }
      }
      if (DEBUG && purePos != -1) {
        std::cout << "Using pure positive column " << purePos << " as pivot\n";
      }

      auto pos = rs.pPlus;
      auto neg = rs.pMinus;
      std::cout << "Preparing to treat row " << targetRow << " with "
          << pos.size () << " positive of which " << pos.basis.size ()
          << " basis and " << neg.size () << " negative columns for a max of "
          << (pos.size () * neg.size ()) << " new columns.\n";
      for (size_t j = 0, je = pos.size (); j < je; ++j) {
        auto jindex = (purePos != -1) ? purePos : pos.keyAt (j);
        if (purePos != -1) {
          j = pos.size ();
        }
        for (size_t k = 0, ke = neg.size (); k < ke; ++k) {
          SparseArray<T> &colj = colsB.getColumn (jindex);
          size_t kindex = neg.keyAt (k);
          SparseArray<T> &colk = colsB.getColumn (kindex);
          T alpha = -colk.get (targetRow);
          T beta = colj.get (targetRow);
          T gcdt = std::gcd (alpha, beta);
          alpha /= gcdt;
          beta /= gcdt;

          if (DEBUG) {
            std::cout << "Computing " << beta << " * " << colk << " + " << alpha
                << " * " << colj << "\n";
          }

          if (purePos != -1 || j == pos.size () - 1) {
            // In-place update
            msut.erase (kindex);  // Erase before modifying colk (mixed-sign)
            auto changed = sumProdInto (beta, colk, alpha, colj);
            normalize (colk);
            if (colk.size () == 0) {
              for (size_t ind = 0; ind < changed.size (); ++ind) {
                size_t key = changed[ind].first;
                rowSigns.setValue (key, kindex, changed[ind].second, false); // Non-basis
              }
              clearColumnWithBasis (kindex, colsB, rowSigns, false); // Non-basis
              if (DEBUG) {
                std::cout << "Cleared empty in-place vector at " << kindex
                    << "\n";
              }
              continue;  // Skip further processing for this column
            }
            if (DEBUG) {
              std::cout << "Built (emplaced at " << kindex << ") " << colk
                  << "\n";
              std::cout << "Changed reported " << changed << "\n";
            }

            if (heur.useQPlusBasis ()) {
              if (colk.isPurePositive ()) {
                // New case: Insert pure positive columns into MSUT
                if (msut.insert (kindex)) {
                  auto [isMinimal, dominated] = hasMinimalSupport (colk, colsB,
                                                                   basisIndices,
                                                                   rowSigns);
                  if (isMinimal) {
                    // Move to basis: clear from non-basis, then set as basis
                    for (size_t ind = 0; ind < changed.size (); ++ind) {
                      size_t key = changed[ind].first;
                      rowSigns.setValue (key, kindex, 0, false); // Clear non-basis
                      rowSigns.setValue (key, kindex, changed[ind].second,
                                         true); // Set basis
                    }
                    for (size_t ind = 0; ind < colk.size (); ++ind) {
                      rowSigns.setValue (colk.keyAt(ind), kindex, 0, false); // Clear non-basis
                      rowSigns.setValue (colk.keyAt(ind), kindex, colk.valueAt(ind), true); // Set basis
                    }
                    for (size_t domIdx : dominated) {
                      msut.erase (domIdx);  // Clear dominated columns from MSUT
                      clearColumnWithBasis (domIdx, colsB, rowSigns, true); // Basis column
                      basisIndices.erase (domIdx);
                      if (DEBUG) {
                        std::cout << "Cleared dominated vector at " << domIdx
                            << "\n";
                      }
                    }
                    basisIndices.insert (kindex);
                    if (DEBUG) {
                      std::cout << "Added minimal in-place vector at " << kindex
                          << ": " << colk << "\n";
                    }
                  } else {
                    msut.erase (kindex);  // Remove non-minimal column from MSUT
                    for (size_t ind = 0; ind < changed.size (); ++ind) {
                      size_t key = changed[ind].first;
                      rowSigns.setValue (key, kindex, changed[ind].second,
                                         false); // Non-basis
                    }
                    clearColumnWithBasis (kindex, colsB, rowSigns, false); // Non-basis
                    filteredCount++;
                    if (DEBUG) {
                      std::cout << "Cleared non-minimal in-place vector at "
                          << kindex << ": " << colk << "\n";
                    }
                  }
                } else {
                  // Duplicate pure positive column detected
                  for (size_t ind = 0; ind < changed.size (); ++ind) {
                    size_t key = changed[ind].first;
                    rowSigns.setValue (key, kindex, 0, false); // Non-basis
                  }
                  clearColumnWithBasis (kindex, colsB, rowSigns, false); // Non-basis
                  if (DEBUG) {
                    std::cout << "Cleared duplicate pure positive vector at "
                        << kindex << ": " << colk << "\n";
                  }
                }
              } else {
                // Existing mixed-sign column logic
                if (msut.insert (kindex)) {
                  for (size_t ind = 0; ind < changed.size (); ++ind) {
                    size_t key = changed[ind].first;
                    rowSigns.setValue (key, kindex, changed[ind].second, false); // Non-basis
                  }
                  if (DEBUG) {
                    std::cout << "Updated unique mixed-sign vector at "
                        << kindex << ": " << colk << "\n";
                  }
                } else {
                  for (size_t ind = 0; ind < changed.size (); ++ind) {
                    size_t key = changed[ind].first;
                    rowSigns.setValue (key, kindex, 0, false); // Non-basis
                  }
                  clearColumnWithBasis (kindex, colsB, rowSigns, false); // Non-basis
                  if (DEBUG) {
                    std::cout << "Cleared duplicate mixed-sign vector at "
                        << kindex << ": " << colk << "\n";
                  }
                }
              }
            } else {
              for (size_t ind = 0; ind < changed.size (); ++ind) {
                size_t key = changed[ind].first;
                rowSigns.setValue (key, kindex, changed[ind].second, false); // Non-basis default
              }
            }
          } else {
            // New column
            auto newCol = sumProd (beta, colk, alpha, colj);
            normalize (newCol);
            if (newCol.size () == 0) {
              if (DEBUG) {
                std::cout << "Skipped empty new vector\n";
              }
              continue;  // Skip appending and further processing
            }
            size_t newColIndex = colsB.getColumnCount ();
            colsB.appendColumn (std::move (newCol));  // Temporarily append
            if (heur.useQPlusBasis ()) {
              if (!colsB.getColumn (newColIndex).isPurePositive ()) {
                if (msut.insert (newColIndex)) {
                  for (size_t ind = 0;
                      ind < colsB.getColumn (newColIndex).size (); ++ind) {
                    rowSigns.setValue (
                        colsB.getColumn (newColIndex).keyAt (ind), newColIndex,
                        colsB.getColumn (newColIndex).valueAt (ind), false); // Non-basis
                  }
                  if (DEBUG) {
                    std::cout << "Built (appended unique mixed-sign at "
                        << newColIndex << ") " << colsB.getColumn (newColIndex)
                        << "\n";
                  }
                } else {
                  colsB.deleteColumn (newColIndex);  // Remove if duplicate
                  if (DEBUG) {
                    std::cout << "Filtered duplicate mixed-sign vector: "
                        << newCol << "\n";
                  }
                }
              } else {
                // New case: Insert pure positive columns into MSUT
                if (msut.insert (newColIndex)) {
                  auto [isMinimal, dominated] = hasMinimalSupport (
                      colsB.getColumn (newColIndex), colsB, basisIndices,
                      rowSigns);
                  if (isMinimal) {
                    // Move to basis: clear from non-basis, then set as basis
                    for (size_t ind = 0;
                        ind < colsB.getColumn (newColIndex).size (); ++ind) {
                      rowSigns.setValue (
                          colsB.getColumn (newColIndex).keyAt (ind),
                          newColIndex, 0, false); // Clear non-basis
                      rowSigns.setValue (
                          colsB.getColumn (newColIndex).keyAt (ind),
                          newColIndex,
                          colsB.getColumn (newColIndex).valueAt (ind), true); // Set basis
                    }
                    for (size_t domIdx : dominated) {
                      msut.erase (domIdx);  // Clear dominated columns from MSUT
                      clearColumnWithBasis (domIdx, colsB, rowSigns, true); // Basis column
                      basisIndices.erase (domIdx);
                      if (DEBUG) {
                        std::cout << "Cleared dominated vector at " << domIdx
                            << "\n";
                      }
                    }
                    basisIndices.insert (newColIndex);
                    if (DEBUG) {
                      std::cout << "Added minimal vector at " << newColIndex
                          << ": " << colsB.getColumn (newColIndex) << "\n";
                    }
                  } else {
                    msut.erase (newColIndex); // Remove non-minimal column from MSUT
                    colsB.deleteColumn (newColIndex);  // Remove if non-minimal
                    filteredCount++;
                    if (DEBUG) {
                      std::cout << "Filtered non-minimal vector: " << newCol
                          << "\n";
                    }
                  }
                } else {
                  // Duplicate pure positive column detected
                  colsB.deleteColumn (newColIndex);  // Remove if duplicate
                  if (DEBUG) {
                    std::cout << "Filtered duplicate pure positive vector: "
                        << newCol << "\n";
                  }
                }
              }
            } else {
              for (size_t ind = 0; ind < colsB.getColumn (newColIndex).size ();
                  ++ind) {
                rowSigns.setValue (colsB.getColumn (newColIndex).keyAt (ind),
                                   newColIndex,
                                   colsB.getColumn (newColIndex).valueAt (ind),
                                   false); // Non-basis
              }
              if (DEBUG) {
                std::cout << "Built (appended at " << newColIndex << ") "
                    << colsB.getColumn (newColIndex) << "\n";
              }
            }
          }
        }
      }
      if (DEBUG) {
        std::cout << "Obtained rowsign " << rowSigns.get (targetRow)
            << " for row " << targetRow << "\n";
      }
      if (rowSigns.get (targetRow).pMinus.size () != 0) {
        std::cout << "Unexpected rowsign " << rowSigns.get (targetRow)
                    << " for row " << targetRow << "\n";
        exit (2);
      }
      posRows.insert (targetRow);

      testPureBasisDominance(colsB, rowSigns, basisIndices, msut, posRows);
    }


    static void testPureBasisDominance (MatrixCol<T> &colsB,
      RowSignsDomination<T> &rowSigns,
      std::unordered_set<size_t> &basisIndices,
      MixedSignsUniqueTable<T> &msut,
      std::unordered_set<size_t> &posRows)
    {
      // now look at basis vectors that use this row => are they only touching pure positive rows ?
      std::vector<size_t> pureBasis;
      for (size_t coli : basisIndices) {
        auto & col = colsB.getColumn (coli);
        for (size_t j=0;j<col.size();j++) {
          if (posRows.find (col.keyAt (j)) == posRows.end ()) {
            break;
          }
          if (j == col.size () - 1) {
            pureBasis.push_back (coli);
          }
        }
      }
      if (! pureBasis.empty()) {
        if (DEBUG) {
          std::cout << "Found pure basis composed of  "<< pureBasis.size () << " vectors\n";
          for (size_t i = 0; i < pureBasis.size (); i++) {
            size_t coli = pureBasis[i];
            auto &col = colsB.getColumn (coli);
            std::cout << "Pure basis vector " << coli << " : " << col << "\n";
          }
        }
        // sort new pure basis by size, increasing
        std::sort (
            pureBasis.begin (), pureBasis.end (),
            [&colsB] (size_t a, size_t b) {
              return colsB.getColumn (a).size () < colsB.getColumn (b).size ();
            });
        std::unordered_set<size_t> toClear;
        // check domination of *any* vector
        for (size_t bindex : pureBasis) {
          const auto &currentCol = colsB.getColumn (bindex);
          auto intersectingCols =
              rowSigns.get (currentCol.keyAt (0)).pPlus.nonBasis;
          // std::cout << "starting inter with " << intersectingCols << "\n";
          for (size_t i = 1; i < currentCol.size (); ++i) {
            size_t row = currentCol.keyAt (i);
            intersectingCols.restrict (
                rowSigns.get (row).pPlus.nonBasis);
            // std::cout << "after refine with " << rowSigns.get (row).pPlus.nonBasis << "current inter " << intersectingCols << "\n";
            if (intersectingCols.size () == 0) {
              break;
            }
          }

          if (DEBUG) {
            std::cout << "Intersection for column " << bindex << " contains "
                << intersectingCols.size () << " columns: ";
            for (size_t j = 0; j < intersectingCols.size (); ++j) {
              std::cout << intersectingCols.keyAt (j) << " ";
            }
            std::cout << "\n";
          }
          if (DEBUG && !intersectingCols.size()==0) {
            std::cout << "Eliminating vectors dominated by pure basis vector " << bindex << ":" << colsB.getColumn(bindex) << " Total "
                << intersectingCols.size () << " columns\n";
          }
          for (size_t k = intersectingCols.size () ; k-->0;) {
            size_t colIdx = intersectingCols.keyAt (k);
            toClear.insert (colIdx);
          }
        }

        if (! toClear.empty()) {
          std::cout << "Eliminating " << toClear.size () << " columns dominated by pure basis vectors\n";
          for (size_t colIdx : toClear) {
            colsB.getColumn (colIdx).clear();
          }
          colsB.dropEmptyColumns ();
          rowSigns = RowSignsDomination<T> (colsB);
          msut.clear ();
          basisIndices.clear ();

          for (size_t i = 0; i < colsB.getColumnCount (); ++i) {
            if (colsB.getColumn (i).isPurePositive ()) {
              basisIndices.insert (i);
            }
            msut.insert (i);
          }
        }

      }


    }

    /**
     * @brief Applies a single row elimination step in Phase 1.
     *
     * Selects a pivot (row, column) using either single-sign row heuristics or general pivot selection, then eliminates
     * the row by updating matC and matB. Tracks the type of elimination (single-sign or general) in counts.
     *
     * @param matC        The incidence matrix (modified).
     * @param matB        The transformation matrix (modified).
     * @param rowSigns    Bookkeeping of row signs (updated).
     * @param counts      Tracks single-sign (first) and general (second) eliminations.
     * @param onlyPositive Enforces semiflow constraints if true.
     * @param heur        Heuristic settings for pivot selection.
     */
    static void applyRowElimination (MatrixCol<T> &matC, MatrixCol<T> &matB,
    RowSigns<T> &rowSigns,
                                     std::pair<size_t, size_t> &counts,
                                     bool onlyPositive,
                                     const EliminationHeuristic &heur)
    {
      // Prioritize single-sign rows for faster elimination when enabled by heuristics.
      // 1) Check Single-Sign rows first:
      if (heur.useSingleSignRow ()) {
        // Single-sign path => pick pivot col from pMinus or pPlus
        // possibly compare matC column sizes if both are size 1
        auto pivot = findSingleSignPivot (matC, rowSigns, heur.getLoopLimit ());
        if (pivot.isSet ()) {
          eliminateRowWithPivot (pivot, matC, matB, rowSigns, onlyPositive);
          counts.first++;
          return;
        }
      }
      // 2) General pivot choice:
      //    (a) pick a column with minimal size in matC
      //    (b) among the rows in that column, pick tRow with the smallest "cost"
      //        (for instance, the smallest absolute cell value, or minimal expansions)
      auto pivot = findBestPivot (matC, rowSigns, onlyPositive,
                                  heur.getLoopLimit ());
      // pivot is set or matC would be zero
      eliminateRowWithPivot (pivot, matC, matB, rowSigns, onlyPositive);
      counts.second++;
    }

  private:

    struct PivotChoice
    {
      ssize_t row;
      ssize_t col;
      PivotChoice (ssize_t r = -1, ssize_t c = -1)
          : row (r), col (c)
      {
      }
      bool isSet () const
      {
        return row != -1;
      }
    };

    static void eliminateRowWithPivot (const PivotChoice pivot,
                                       MatrixCol<T> &matC, MatrixCol<T> &matB,
                                       RowSigns<T> &rowSigns,
                                       bool onlyPositive)
    {
      assert(pivot.isSet ());
      const auto &tRow = pivot.row;
      const auto &tCol = pivot.col;
      const auto &rowsign = rowSigns.get (tRow);

      SparseBoolArray toVisit = SparseBoolArray::unionOperation (rowsign.pMinus,
                                                                 rowsign.pPlus);

      if (toVisit.size () == 0) {
        // nothing to do
        return;
      }
      if (onlyPositive
          && (rowsign.pMinus.size () == 0 || rowsign.pPlus.size () == 0)) {
        // [1.1] if there exists a row h in C such that the sets P+ = {j | c_hj > 0},
        // P- = {j | c_hj < 0} satisfy P+ == {} or P- == {} and not (P+ == {} and P- == {})
        // that means it exists a row that all components are positive respectively negative
        // [1.1.a] delete from the extended matrix all the columns of index j \in P+ \cup P-
        for (size_t i = 0, ie = toVisit.size (); i < ie; i++) {
          clearColumn (toVisit.keyAt (i), matC, matB, rowSigns);
        }
        return;
      }

      T cHk = matC.get (tRow, tCol);
      T bbeta = std::abs (cHk);

      if (DEBUG) {
        std::cout << "tCol : " << tCol << " tRow " << tRow << std::endl;
        std::cout << "rowsign : " << rowsign << std::endl;
      }

      for (size_t i = 0; i < toVisit.size (); i++) {
        size_t j = toVisit.keyAt (i);
        if (j == (size_t) tCol) {
          continue;
        }

        SparseArray<T> &colj = matC.getColumn (j);
        T cHj = colj.get (tRow);
        if (cHj != 0) {
          // Compute alpha and beta
          T alpha =
              ((signum (cHj) * signum (cHk)) < 0) ?
                  std::abs (cHj) : -std::abs (cHj);
          if (alpha == 0 && bbeta == 1) {
            // No operation needed
            continue;
          }
          T gcdt = std::gcd (alpha, bbeta);
          alpha /= gcdt;
          T beta = bbeta / gcdt;

          if (DEBUG) {
            std::cout << "Eliminating with pivot col=" << tCol << " row="
                << tRow << " combining j=" << j << std::endl;
          }

          // Update matC
          auto changed = sumProdInto (beta, colj, alpha, matC.getColumn (tCol));
          // Update the row-sign bookkeeping
          for (size_t ind = 0, inde = changed.size (); ind < inde; ind++) {
            size_t key = changed[ind].first;
            rowSigns.setValue (key, j, changed[ind].second);
          }

          // Update matB
          SparseArray<T> &coljb = matB.getColumn (j);
          sumProdIntoNoChange (beta, coljb, alpha, matB.getColumn (tCol));

          // Optionally perform GCD reduction on column j
          T gcdm = gcd (matC.getColumn (j));
          if (gcdm != 1) {
            T gcdb = gcd (matB.getColumn (j));
            if (gcdb != 1) {
              T ggcd = std::gcd (gcdm, gcdb);
              if (ggcd != 1) {
                matC.getColumn (j).scalarDiv (ggcd);
                matB.getColumn (j).scalarDiv (ggcd);
              }
            }
          }
        }
      }
      // Finally clear the pivot column tCol
      clearColumn (tCol, matC, matB, rowSigns);
    }

    static PivotChoice findSingleSignPivot (const MatrixCol<T> &matC,
                                            const RowSigns<T> &rowSigns,
                                            size_t loopLimit)
    {

      PivotChoice pivot;
      ssize_t candidateRow = rowSigns.findSingleSignRow (loopLimit);
      if (candidateRow != -1) {
        // Single-sign path => pick pivot col from pMinus or pPlus
        // possibly compare matC column sizes if both are size 1
        // Get the candidate row data freshly.
        const auto &rowData = rowSigns.get (candidateRow);
        // In our construction, exactly one of pPlus or pMinus must have size 1.
        assert(rowData.pPlus.size() == 1 || rowData.pMinus.size() == 1);

        // Determine which set is unique.
        // If pPlus is unique then tCol is the unique key and the complementary set is pMinus;
        // otherwise, tCol is from pMinus and the complementary set is pPlus.
        bool isNeg = (rowData.pMinus.size () == 1);
        ssize_t tCol =
            isNeg ? rowData.pMinus.keyAt (0) : rowData.pPlus.keyAt (0);

        if (rowData.pPlus.size () == 1 && rowData.pMinus.size () == 1) {
          // std::cout << "Examine col j=" << j << " size=" << matC.getColumn (j).size () <<" with col=" << tCol << " size " << matC.getColumn (tCol).size () << std::endl;
          // we can actually choose which one to get rid of
          if (matC.getColumn (rowData.pMinus.keyAt (0)).size ()
              > matC.getColumn (rowData.pPlus.keyAt (0)).size ()) {
            tCol = rowData.pPlus.keyAt (0);
            isNeg = !isNeg;
            // std::cout << "SWAPPED" << std::endl;
          }
        }

        return PivotChoice
          { candidateRow, tCol };
      } else {
        return PivotChoice ();
      }
    }

    template<typename RS>
      static PivotChoice findBestPivot (const MatrixCol<T> &matC,
                                        const RS &rowSigns, bool onlyPositive,
                                        size_t loopLimit)
      {
        // We'll assume there's at least one non-empty column; otherwise we can't pivot.
        // If you want to handle an all-empty matrix, you can do so by checking further.

        size_t bestColSize = std::numeric_limits<size_t>::max ();
        std::vector<size_t> candidateCols;
        candidateCols.reserve (matC.getColumnCount ());

        // 1) Find the columns that have the minimal nonzero size in one pass.
        for (size_t c = 0; c < matC.getColumnCount (); c++) {
          size_t s = matC.getColumn (c).size ();
          if (s == 0) {
            continue; // skip empty columns
          }
          if (s < bestColSize) {
            bestColSize = s;
            candidateCols.clear ();
            candidateCols.push_back (c);
          } else if (s == bestColSize) {
            candidateCols.push_back (c);
          }
          if (candidateCols.size () >= loopLimit) {
            break; // early exit if we have enough candidates
          }
        }
        // Here we have a list of columns whose size is 'bestColSize'.

        // 2) Among those columns, find the best (row,col) pivot by row heuristics:
        //    - smallest rowSize = pPlus.size() + pMinus.size()
        //    - tie-break on absolute value
        //    - early break if rowSize <= 2 and absVal == 1
        ssize_t bestRow = 0;
        ssize_t bestCol = candidateCols[0];       // default
        size_t bestRowSz = std::numeric_limits<size_t>::max ();
        T bestAbsVal = std::numeric_limits<T>::max ();
        bool foundPivot = false;

        // For each candidate column, iterate over its nonzero entries
        for (size_t c : candidateCols) {
          const SparseArray<T> &colData = matC.getColumn (c);
          for (size_t i = 0, ie = colData.size (); i < ie; i++) {
            size_t r = colData.keyAt (i);
            T val = colData.valueAt (i);
            // We assume val != 0 in a properly stored SparseArray.

            // Row size in rowSigns
            const auto &rs = rowSigns.get (r);

            if (onlyPositive
                && (rs.pMinus.size () == 0 || rs.pPlus.size () == 0)) {
              // We can't pivot on a row that has all positive or all negative entries.
              return PivotChoice (r, c);
            }
            size_t rowSz = rs.pPlus.size () + rs.pMinus.size ();
            T absVal = (val < 0) ? -val : val;

            // Compare to see if it's better:
            //   - pick the pivot that has the smaller rowSz
            //   - tie-break on smaller absVal
            if (!foundPivot || (rowSz < bestRowSz)
                || (rowSz == bestRowSz && absVal < bestAbsVal)) {
              foundPivot = true;
              bestRow = r;
              bestCol = c;
              bestRowSz = rowSz;
              bestAbsVal = absVal;

              // Early‐exit if (rowSz <= 2) and (absVal == 1).
              if (rowSz <= 2 && absVal == 1) {
                // Immediately return this pivot.
                return PivotChoice
                  { bestRow, bestCol };
              }
            }
          }
        }

        // If we found no pivot at all, that means all columns might have size 0
        // or we never had nonzero entries. But in practice we do handle above.
        // So we finalize the best we found:
        return foundPivot ? PivotChoice
          { bestRow, bestCol } :
                            PivotChoice ();  // an unset pivot
      }

    template<typename RS>
      static ssize_t findBestFMERow (const MatrixCol<T> &matB,
                                     const RS &rowSigns, size_t loopLimit)
      {
        ssize_t minRow = -1;
        ssize_t minRowWeight = std::numeric_limits<ssize_t>::max ();
        ssize_t minRefinedRowWeight = std::numeric_limits<ssize_t>::max ();
        size_t seen = 0;
        for (const auto &rs : rowSigns) {
          size_t pps = rs.pPlus.size ();
          size_t ppm = rs.pMinus.size ();

          if (seen++ > loopLimit && minRow != -1) {
            break;
          }

          if (pps == 0)
          // a pure negative row => we need to clear all related columns
          return rs.row;

          // number of columns in result
          ssize_t weight = pps * ppm - pps - ppm;

          if (minRow == -1 || minRowWeight > weight) {
            ssize_t refinedweight = 0;
            for (size_t i = 0, ie = rs.pPlus.size (); i < ie; i++) {
              refinedweight += matB.getColumn (rs.pPlus.keyAt (i)).size ();
            }
            for (size_t i = 0, ie = rs.pMinus.size (); i < ie; i++) {
              refinedweight += matB.getColumn (rs.pMinus.keyAt (i)).size ();
            }
            if (minRow == -1 || minRefinedRowWeight > refinedweight) {
              minRow = rs.row;
              minRowWeight = weight;
              minRefinedRowWeight = refinedweight;
            }
          }
        }

        if (DEBUG) {
          if (minRow != -1) {
            std::cout << "Selected row " << minRow << " with weight "
                << minRowWeight << " and refined weight " << minRefinedRowWeight
                << std::endl;
          } else {
            std::cout << "No row selected from " << rowSigns << " with matB= "
                << matB << std::endl;
          }
        }

        return minRow;
      }

    template<typename RS>
      static std::pair<bool, std::vector<size_t>> hasMinimalSupport (
          const SparseArray<T> &newCol, const MatrixCol<T> &basis,
          const std::unordered_set<size_t> &basisIndices, const RS &rowSigns)
      {
        std::vector<size_t> counts (basis.getColumnCount (), 0); // Pre-allocate
        size_t newColSize = newCol.size ();

        // Check if newCol is dominated
        for (size_t i = 0; i < newColSize; ++i) {
          size_t row = newCol.keyAt (i);
          const auto &rs = rowSigns.get (row);
          for (size_t j = 0, je = rs.pPlus.basis.size (); j < je; ++j) {
            size_t colIdx = rs.pPlus.basis.keyAt (j);

            if (++counts[colIdx] == basis.getColumn (colIdx).size ()
                && basisIndices.count (colIdx) > 0) {
              if (DEBUG) {
                std::cout << "Vector " << newCol
                    << " dominated by basis column " << colIdx << " ("
                    << basis.getColumn (colIdx).size () << ")\n";
              }
              return {false, {}};
            }
          }
        }

        // Check if newCol dominates any basis vectors
        std::vector<size_t> dominated;
        for (size_t colIdx : basisIndices) {
          if (counts[colIdx] == newColSize
              && basis.getColumn (colIdx).size () > newColSize) {
            if (DEBUG) {
              std::cout << "Vector " << newCol << " dominates basis column "
                  << colIdx << " (" << basis.getColumn (colIdx).size ()
                  << ")\n";
            }
            dominated.push_back (colIdx);
          }
        }

        return {true, dominated};
      }

    static std::pair<std::vector<Permutation>, MatrixCol<T>> factorizeBasis (
        MatrixCol<T> &colsB)
    {
      auto startTime = std::chrono::high_resolution_clock::now ();

      // Step 0: Normalize input
      std::cout << "Normalizing flow basis matrix with "
          << colsB.getColumnCount () << " columns\n";
      colsB.normalizeAndReduce (true);
      std::cout << "Produced matrix with " << colsB.getColumnCount ()
          << " columns\n";

      // Step 1: Sort by size (descending, largest first)
      colsB.sortByColumnSize (true);

      // Step 2: Build row index, basic version no indexes, no positive treatment
      RowSigns<T> rowSigns (colsB, false, false);

      // Step 3: Factorize
      std::vector<Permutation> permutations;

      std::cout << "Starting factorization with " << colsB.getColumnCount ()
          << " columns\n";
      prefilterAndFactorize (colsB, rowSigns, permutations);

      // Step 4: Finalize reduced basis
      if (permutations.size () > 0) colsB.normalizeAndReduce (true);

      // Step 5: Report
      auto endTime = std::chrono::high_resolution_clock::now ();
      auto duration = std::chrono::duration_cast < std::chrono::milliseconds
          > (endTime - startTime);
      std::cout << "Factorization complete. Permutations: "
          << permutations.size () << ", Reduced basis size: "
          << colsB.getColumnCount () << ". Time: " << duration.count ()
          << " ms\n";

      fusePermutations(permutations,colsB);

      return {permutations, colsB};
    }


    static void fusePermutations (std::vector<Permutation> &permutations,
                                  MatrixCol<T> &colsB)
    {
      std::unordered_map<size_t, size_t> indexMap;
      for (size_t i = 0; i < permutations.size (); i++) {
        indexMap[permutations[i].index] = i;
      }

      bool needRewrite = false;
      //iterate permutations in reverse order
      for (size_t i = permutations.size (); i-- > 0;) {
        auto &p = permutations[i];
        // iterate parts
        bool redo = true;
        while (redo) {
          redo = false;
          for (size_t j = 0; j < p.elements.size() ; j++) {
            auto &e = p.elements[j];
            // check condition : at most one expansion per term
            int nbterms =0;
            for (size_t k = 0; k < e.size (); k++) {
              auto key = e.keyAt (k);
              auto it = indexMap.find (key);
              if (it != indexMap.end ()) {
                nbterms++;
                if (nbterms > 1) {
                  break;
                }
              }
            }
            if (nbterms > 1) {
              // skip these expansions
              continue;
            }
            // for each key, is it in our set ?
            for (size_t k = 0; k < e.size (); k++) {
              auto key = e.keyAt (k);
              auto it = indexMap.find (key);
              if (it != indexMap.end ()) {
                auto copy = e;
                p.elements.erase (p.elements.begin () + j);
                T coeff = copy.valueAt (k);
                copy.put (key, 0);
                for (auto &rep : permutations[it->second].elements) {
                    p.elements.emplace_back(sumProd((T)1, copy, coeff, rep));
                }
                redo = true;
                needRewrite = true;
                break;
              }
            }
            if (redo)
              break;
          }
        }
      }


      if (needRewrite) {
        size_t pindex = colsB.getRowCount() - permutations.size();
        std::vector<Permutation> newperms;
        auto trans = colsB.transpose ();
        std::vector<size_t> todel;
        for (auto & p : permutations) {
          if (trans.getColumn(p.index).size () == 0) {
            todel.push_back(p.index);
          } else {
            p.index = pindex++;
            newperms.push_back (p);
          }
        }
        permutations = newperms;
        std::sort (todel.begin (), todel.end (), std::greater<size_t>());
        for (auto &d : todel) {
          trans.deleteColumn (d);
        }
        colsB = trans.transpose();
      }
    }

    static void prefilterAndFactorize (MatrixCol<T> &colsB,
    RowSigns<T> &rowSigns,
                                       std::vector<Permutation> &permutations)
    {
      std::vector<bool> isNonFactorizable (colsB.getColumnCount (), false);

      for (size_t fIdx = colsB.getColumnCount (); fIdx-- > 0;) {
        if (isNonFactorizable[fIdx]) continue;
        const SparseArray<T> &F = colsB.getColumn (fIdx);

        // Skip fully positive vectors
        if (F.size () == 0 || F.isPurePositive ()) continue;

        std::unordered_set<size_t> targets;
        // Trivial equality check : size=2 and not pure positive and weights are +-1
        if (F.size () == 2 && std::abs (F.valueAt (0)) == 1
            && std::abs (F.valueAt (1)) == 1) {
          if (DEBUG) {
            std::cout << "Trivially factorizable flow " << F << " at index "
                << fIdx << "\n";
          }
          // these ones are ALWAYS factorizable
          for (size_t j = 0; j < 2; ++j) {
            const auto &rs = rowSigns.get (F.keyAt (j));
            for (size_t i = 0; i < rs.pPlus.size (); ++i) {
              targets.insert (rs.pPlus.keyAt (i));
            }
            for (size_t i = 0; i < rs.pMinus.size (); ++i) {
              targets.insert (rs.pMinus.keyAt (i));
            }
          }
        } else {
          targets = prefilterFactorizability (fIdx, colsB, rowSigns,
                                              isNonFactorizable);
        }
        if (targets.empty () || isNonFactorizable[fIdx]) continue;

        if (DEBUG) {
          std::cout << "Prefiltered to obtain " << targets.size ()
              << " possible factorizable columns with F " << F << "at " << fIdx
              << "\n";
          std::cout << "Targets :\n";
          for (auto &t : targets) {
            std::cout << t << " ";
          }
          std::cout << std::endl;
        }

        // at this stage, F is at least possibly factorizable
        // and we have a list of targets to check for existence of coefficients k+,k-

        // to store the k+,k- we find
        std::map<size_t, std::pair<Rational<T>, Rational<T>>> coefficients;
        for (const auto &target : targets) {
          auto [success, kpkm] = SparseArray<T>::isFactorizable (
              F, colsB.getColumn (target));
          if (success) {
            coefficients[target] = kpkm;
          } else {
            isNonFactorizable[target] = true;
            isNonFactorizable[fIdx] = true;
          }
        }
        // coefficients didn't line up
        if (isNonFactorizable[fIdx]) continue;

        if (DEBUG) {
          std::cout << "Before rewriting :" << colsB << std::endl;
          std::cout << "Targets :";
          for (auto &t : targets) {
            std::cout << t << " ";
          }

          std::cout << "\nCoefficients :\n";
          for (const auto& [j, kpkm] : coefficients) {
            std::cout << j << " : " << kpkm.first << " " << kpkm.second
                << std::endl;
          }
        }
        // Factorize
        permutations.push_back (
            rewriteBasis (fIdx, coefficients, colsB, rowSigns));

        if (DEBUG) {
          std::cout << "After rewriting using "
              << permutations[permutations.size () - 1] << " obtained :\n"
              << colsB << std::endl;
        }

      }
    }

    /**
     * @brief Rewrites the basis by factorizing flow F at fIdx, introducing a new representative.
     * @param fIdx Index of flow F in colsB to factorize out.
     * @param coefficients Map of column indices to (k+, k-) pairs for each affected flow F', computed earlier.
     * @param colsB Basis of flows as sparse column vectors, modified in-place to reflect factorization.
     * @param rowSigns Sign index with pPlus (positive) and pMinus (negative) per row, updated with new rep.
     * @return Permutation object with new row index and F+/F- as positive vectors.
     *
     * Given a factorizable F, for each F' in colsB F' + k+ F+ + k- F- = A holds:
     * - A is F' outside F’s support (A[i] = F'[i] if i not in |F|, else 0).
     * - We introduce rep_i = F+ - F- and rewrite F' = A - (k+ + k-) rep_i.
     * - F itself becomes rep_i - rep_i = 0, nullified and cleared.
     * Updates colsB and rowSigns accordingly, adding rep_i as a new row and
     * returning the new permutation induced by F.
     */
    static Permutation rewriteBasis (
        ssize_t fIdx, std::map<size_t, std::pair<Rational<T>, Rational<T>>> &coefficients,
        MatrixCol<T> &colsB, RowSigns<T> &rowSigns)
    {
      // F is the factorizable flow to rewrite out of the basis.
      const SparseArray<T> &F = colsB.getColumn (fIdx);

      // NB: all conditions are met, "coefficients"'s keys are the vector
      // indexes we need to rewrite, and each comes with its k+/k- coefficients.
      if (DEBUG) {
        std::cout << "Found factorizable flow " << F << " at index " << fIdx
            << std::endl;
      }

      // Collect the F+ and F- components to prepare the permutation object we return.
      // We need to copy them before clearing all those rows in the whole matrix
      SparseArray<T> Fplus, Fminus;
      // rows to delete : support of F
      std::vector<size_t> toDel;
      for (size_t i = 0, ie = F.size (); i < ie; ++i) {
        size_t row = F.keyAt (i);
        T val = F.valueAt (i);
        if (val > 0) Fplus.append (row, val);
        else Fminus.append (row, -val);
        toDel.push_back (row);
      }

      // Create permutation with F+ and F- as its terms
      // and as index a newly created row in the basis
      size_t newIdx = colsB.appendRow ();
      Permutation perm (newIdx,
        { Fplus, Fminus });

      // Rewrite basis
      // 1. Clear all rows in F

      // sort todel greatest to smallest index :
      // deleting from end (higher indexes) is more
      // efficient on SparseArray, less shifting.
      std::sort (toDel.begin (), toDel.end (),
                 [] (auto &a, auto &b) {return a > b;});
      for (size_t i = 0, ie = colsB.getColumnCount (); i < ie; ++i) {
        auto &col = colsB.getColumn (i);
        for (const auto &row : toDel) {
          col.put (row, 0);
        }
      }
      // Also cleanup RowSigns, which must stay synchronized with colsB
      // for the next call
      for (auto row : toDel) {
        rowSigns.clearRow (row);
      }

      // Step 2: Rewrite each F' = lcm(k+,k-) * (A - (k+ + k-) rep_i).
      // A is F' outside |F|, already handled by clearing |F|.
      for (const auto& [j, kpkm] : coefficients) {
        SparseArray<T> &Fprime = colsB.getColumn (j);

        // Total coefficient for rep_i.
        T lcm = std::lcm (kpkm.first.denominator, kpkm.second.denominator);
        // rescale k+ and k- to integers
        T k_plus = kpkm.first.numerator * lcm / kpkm.first.denominator;
        T k_minus = kpkm.second.numerator * lcm / kpkm.second.denominator;
        T k = -k_plus + k_minus;
        // No change if k+ + k- = 0 ; This happens for instance for F itself.
        if (k == 0) continue;

        if (DEBUG) {
          std::cout << "Rewriting flow " << Fprime << " with k+ " << k_plus
              << " and k-" << k_minus  << " and lcm " << lcm << " using rep " << newIdx
              << std::endl;
        }

        // multiply by lcm
        if (lcm != 1) {
          Fprime.scalarProduct(lcm);
        }
        // Add -(k+ + k-) at rep_i’s row, completing F' = A - (k+ + k-) rep_i.
        Fprime.append (newIdx, -k);
        normalizeWithSign(Fprime);
        // keep row index in sync
        rowSigns.setValue (newIdx, j, -k);
      }

      // check nullified F
      if (colsB.getColumn (fIdx).size () != 0) {
        std::cerr << "Error: Factorization failed, F not nullified\n"
            << colsB.getColumn (fIdx) << std::endl;
        assert(false);
      }

      return perm;
    }

    /**
     * @brief Prefilters columns that could be factorized with flow F at fIdx.
     * @param fIdx Index of flow F in colsB to test against others.
     * @param colsB Basis of flows as sparse column vectors (modified in-place elsewhere).
     * @param rowSigns Sign index with pPlus (positive) and pMinus (negative) per row.
     * @param isNonFactorizable Boolean vector, set to true for non-factorizable columns detected during the prefiltering.
     * @return A set of column indices that might factorize with F, empty if F is not factorizable.
     * Coefficient values are not checked, so these columns still need to pass the true SparseArray::isFactorizable test.
     * Also updates the tags in isNonFactorizable (sets some to true).
     *
     * Checks sign consistency and full overlap of other flows F' with F+ and F-.
     * Uses pos/neg counters to track polarity and detect conflicts.
     */
    static std::unordered_set<size_t> prefilterFactorizability (
        size_t fIdx, MatrixCol<T> &colsB,
        RowSigns<T> &rowSigns,
        std::vector<bool> &isNonFactorizable)
    {
      // F is the target candidate flow to factorize, extract from colsB.
      const SparseArray<T> &F = colsB.getColumn (fIdx);

      // pos, neg: count intersections of all flows F' with the F+ (resp. F-) component of F.
      // Positive value means F' entries are positive, negative means negative.
      // preallocated and zeroed out to number of flows in colsB.
      std::vector<int> pos (colsB.getColumnCount (), 0);
      std::vector<int> neg (colsB.getColumnCount (), 0);

      // to compute the sizes of |F+| and |F-| as we go
      size_t fplusSize = 0, fminusSize = 0;

      // Count overlaps of F with other flows F' in colsB
      for (size_t i = 0, ie = F.size (); i < ie; ++i) {
        size_t row = F.keyAt (i);
        bool isPos = F.valueAt (i) > 0;

        // update our size counts for F+/F-
        if (isPos) {
          fplusSize++;
        } else {
          fminusSize++;
        }

        // depending on isPos, we inc/decrement the pos or neg counters
        auto &pn = isPos ? pos : neg;

        // the current row's sign information
        // efficient per row access
        const auto &rs = rowSigns.get (row);

        // iterate over both parts of rowsign
        // so first rs.pPlus then rs.pMinus
        for (size_t j = 0, je = rs.pPlus.size (); j < je; ++j) {
          // flow F'
          size_t colIdx = rs.pPlus.keyAt (j);
          // before incrementing, check that we don't have a conflict
          if (pn[colIdx] < 0) {
            // Conflict: F' was negative before, now positive.
            // Example: F = p0 + p1 - p2, F' = p0 - p1 +...
            // so signs of p0 and p1 in F are inconsistent with F'+

            // Both vectors are now tagged as non-factorizable :
            // its a global property wrt to the whole basis that isNonFactorizable tracks.
            isNonFactorizable[fIdx] = true;
            isNonFactorizable[colIdx] = true;

            // no break, we want to tag more columns F' if we can
            // using the current (small at start) F
          } else {
            // update the counter, we intersect F+ => pos
            // and we are in rs.pPlus so we increment
            pn[colIdx]++;
          }
        }
        // symmetric/dual case for pMinus
        for (size_t j = 0, je = rs.pMinus.size (); j < je; ++j) {
          size_t colIdx = rs.pMinus.keyAt (j);
          if (pn[colIdx] > 0) {
            isNonFactorizable[colIdx] = isNonFactorizable[fIdx] = true;
          } else {
            pn[colIdx]--;
          }
        }
      }

      // no break or return, we still want to tag more columns F' if we can
      // we've only tagged sign conflicts so far

      // Check partial overlaps
      // it's forbidden to have a partial overlap with |F+| or |F-|
      // you either contain it whole or not at all.
      std::unordered_set<size_t> targets;
      for (size_t i = 0, ie = colsB.getColumnCount (); i < ie; ++i) {
        // Check if the target is a candidate for factorization
        size_t p = std::abs (pos[i]), n = std::abs (neg[i]);
        // no intersection
        if (p == 0 && n == 0) continue;
        // bad intersection
        if ((p != 0 && p != fplusSize) || (n != 0 && n != fminusSize)) {
          // this is a partial overlap
          // due to the way we don't increment counters in sign conflicts
          // it will also (re) tag them
          isNonFactorizable[i] = isNonFactorizable[fIdx] = true;
        } else {
          // F' matches our "full containment or no overlap" requirement
          // with both F+ and F-
          targets.insert (i);
        }
      }

      if (isNonFactorizable[fIdx]) {
        if (DEBUG) {
          std::cout << "Non-factorizable flow " << F << " at index " << fIdx
              << std::endl;
        }
        return {};
      }

      // a side effect of the test is that F's index fIdx itself is selected by this filter
      // this is by design, it will factorize itself out to 0 if we actually
      // go through with factorizing F
      return targets;
    }

  public:

    static void clearColumn (size_t tCol, MatrixCol<T> &matC,
                             MatrixCol<T> &matB,
                             RowSigns<T> &rowSigns)
    {
      // delete from the extended matrix the column of index k
      SparseArray<T> &colk = matC.getColumn (tCol);
      for (size_t i = 0, ie = colk.size (); i < ie; i++) {
        rowSigns.setValue (colk.keyAt (i), tCol, 0);
      }
      colk.clear ();
      matB.getColumn (tCol).clear ();
    }

    template<typename RS>
      static void clearColumn (size_t tCol, MatrixCol<T> &matB, RS &rowSigns)
      {
        // delete from the extended matrix the column of index k
        SparseArray<T> &colk = matB.getColumn (tCol);
        for (size_t i = 0, ie = colk.size (); i < ie; i++) {
          rowSigns.setValue (colk.keyAt (i), tCol, 0);
        }
        colk.clear ();
      }

    template<typename RS>
      static void clearColumnWithBasis (size_t tCol, MatrixCol<T> &matB,
                                        RS &rowSigns, bool inBasis)
      {
        // Delete from the extended matrix the column of index tCol
        SparseArray<T> &colk = matB.getColumn (tCol);
        for (size_t i = 0, ie = colk.size (); i < ie; i++) {
          rowSigns.setValue (colk.keyAt (i), tCol, 0, inBasis);
        }
        colk.clear ();
      }

  };

}

#endif /* INVARIANTCALCULATOR_H_ */
