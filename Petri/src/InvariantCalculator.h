#ifndef INVARIANTCALCULATOR_H_
#define INVARIANTCALCULATOR_H_

/**
 * A calculator for invariants and testing if a net is covered by invariants.
 * Provides two differient algorithms for calculating the invariants. The first
 * algorithm is descripted in http://de.scribd.com/doc/49919842/Pn-ESTII (slide
 * 88) and the other is also based an the farkas algorithm and is descripted in
 * http://pipe2.sourceforge.net/documents/PIPE-Report.pdf (page 19) which is
 * based on the paper of D'Anna and Trigila "Concurrent system analysis using
 * Petri nets – an optimised algorithm for finding net invariants", Mario D'Anna
 * and Sebastiano Trigila, Computer Communications vol 11, no. 4 august 1988.
 *
 * @author Dennis-Michael Borde, Manuel Gieseking , Adapted to ITS-tools by Yann
 *         Thierry-Mieg, 2017.
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
#include "InvariantsTrivial.h"
#include "Heuristic.h"

namespace petri
{

template<typename T>
  class InvariantCalculator
  {

    static inline const bool DEBUG = false;

  public:
    /**
     * Enumeration for choosing which algorithm should be used.
     */
    enum class InvariantAlgorithm
    {
      // Farkas,
      PIPE
    };

  private:
    /**
     * Hidden constructor
     */
    InvariantCalculator ()
    {
    }

    /**
     * Calculates the invariants with the algorithm based on
     * http://pipe2.sourceforge.net/documents/PIPE-Report.pdf (page 19).
     *
     * @param mat          - the matrix to calculate the invariants from.
     * @param onlyPositive whether we just stop at Flows or go for Semi-Flows
     * @param pnames       variable names
     * @return a generator set of the invariants.
     */
  public:
    static std::unordered_set<SparseArray<T>> calcInvariantsPIPE (
        MatrixCol<T> mat, bool onlyPositive, const EliminationHeuristic &heur =
            EliminationHeuristic ())
    {
      if (mat.getColumnCount () == 0 || mat.getRowCount () == 0) {
        return std::unordered_set<SparseArray<T>> ();
      }
      MatrixCol<T> tmat = mat.transpose ();
      std::unordered_set<SparseArray<T>> normed = std::unordered_set<
          SparseArray<T>> ();
      for (size_t i = 0; i < tmat.getColumnCount (); i++) {
        SparseArray<T> &norm = tmat.getColumn (i);
        normalize (norm);
        normed.insert (norm);
      }
      if (normed.size () < tmat.getColumnCount ()) {
        std::cout << "Normalized transition count is " << normed.size ()
            << " out of " << tmat.getColumnCount () << " initially."
            << std::endl;
      }
      MatrixCol<T> matnorm (tmat.getRowCount (), 0);
      for (const SparseArray<T> &col : normed) {
        matnorm.appendColumn (col);
      }

      MatrixCol<T> matB = phase1PIPE (matnorm.transpose (), onlyPositive, heur);

//		const MatrixCol<T> matB = phase1PIPE(new MatrixCol<T>(mat));
      // We want to work with columns in this part of the algorithm
      // We add and remove columns all day => we want to switch to a column based
      // representation
      // order of rows is really irrelevant + columns which are identical up to
      // scaling factor are useless
      // let's use a set of columns.
      std::unordered_set<SparseArray<T>> semiFlows (2 * matB.getColumnCount ());
      for (size_t i = 0; i < matB.getColumnCount (); i++) {
        SparseArray<T> &col = matB.getColumn (i);
        if (col.size () != 0) {
          normalizeWithSign (col);
          semiFlows.insert (col);
        }
      }

      if (!onlyPositive) {
        return semiFlows;
      }

      MatrixCol<T> colsB (tmat.getRowCount (), 0);
      for (const SparseArray<T> &cb : semiFlows) {
        colsB.appendColumn (cb);
      }
      semiFlows.clear ();
      // phase 2
      std::cout << "// Phase 2 : computing semi flows from basis of "
          << colsB.getColumnCount () << " invariants " << std::endl;

      phase2Pipe (colsB, semiFlows, heur);

      return semiFlows;
    }

  private:

    static MatrixCol<T> phase1PIPE (MatrixCol<T> matC, bool onlyPositive,
                                    const EliminationHeuristic &heur)
    {
      // Build the initial transformation matrix.
      MatrixCol<T> matB = MatrixCol<T>::identity (matC.getColumnCount (),
                                                  matC.getColumnCount ());

      // Vector to hold trivial invariants.
      std::vector<SparseArray<T>> trivialInv;
      // Remove trivial invariants (empty columns in matC) early.
      cullConstantColumns (matC, matB, trivialInv);
      // Remove duplicate columns
      cullDuplicateColumns (matC, matB, trivialInv);
      std::cout << "// Phase 1: matrix " << matC.getRowCount () << " rows "
          << matC.getColumnCount () << " cols " << matC.getEntryCount ()
          << " entries" << std::endl;
      RowSigns rowSigns (matC, heur.useSingleSignRow ());

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

      std::pair<size_t, size_t> counts (0, 0);
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

    static void phase2Pipe (MatrixCol<T> &colsB,
                            std::unordered_set<SparseArray<T>> &semiFlows,
                            const EliminationHeuristic &heur)
    {
      RowSigns<T> rowSigns (colsB, true, true);
      // step 1 : clear pure negative rows.
      {
        std::vector<size_t> tokill;
        for (const auto &rs : rowSigns) {
          if (rs.pPlus.size () == 0) {
            for (size_t i = 0, ie = rs.pMinus.size (); i < ie; i++) {
              colsB.getColumn (rs.pMinus.keyAt (i)).clear ();
            }
            tokill.push_back (rs.row);
          }
        }
        for (const auto &row : tokill) {
          rowSigns.clearRow (row);
        }
      }

      while (true) {
        ssize_t tRow = rowSigns.findSingleSignRow (heur.getLoopLimit ());

        if (tRow == -1) {
          // look for a "small" row
          tRow = findBestFMERow (colsB, rowSigns, heur.getLoopLimit ());
        }

        if (tRow == -1) {
          break;
        }

        eliminateRowFME (tRow, colsB, rowSigns);
      }

      // all coefficients are positive !
      for (size_t i = 0; i < colsB.getColumnCount (); i++) {
        auto &col = colsB.getColumn (i);
        if (col.size () > 0) {
          normalizeWithSign (col);
          semiFlows.insert (col);
        }
      }

      // all done
    }

    static void eliminateRowFME (ssize_t targetRow, MatrixCol<T> &colsB,
    RowSigns<T> &rowSigns)
    {
      const auto &rs = rowSigns.get (targetRow);

      if (DEBUG) {
        std::cout << "Eliminating row " << targetRow << " with "
            << rs.pPlus.size () << " plus and " << rs.pMinus.size ()
            << " minus columns" << std::endl;
        std::cout << rs << std::endl;
      }

      if (rs.pPlus.size () == 0) {
        // this whole row sucks, all columns touching it suck too.
        for (size_t i = 0, ie = rs.pMinus.size (); i < ie; i++) {
          colsB.getColumn (rs.pMinus.keyAt (i)).clear ();
        }
        rowSigns.clearRow (targetRow);
        std::cout << "Cleared row " << targetRow << std::endl;
        return;
      }
      ssize_t purePos = -1;
      if (rs.pPlus.size () > 1) {
        // we might be lucky to have (small) pure positive column in P+
        // We can use it as single pivot to eliminate P- in the row
        // it's guaranteed at least to not introduce new negative values in the columns we touch
        for (size_t i = 0, ie = rs.pPlus.size (); i < ie; i++) {
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
        std::cout << "Using pure positive column " << purePos << " as pivot"
            << std::endl;
      }

      for (size_t j = 0, je = rs.pPlus.size (); j < je; j++) {
        auto jindex = rs.pPlus.keyAt (j);
        if (purePos != -1) {
          jindex = purePos;
          j = je;
        }
        for (size_t k = 0, ke = rs.pMinus.size (); k < ke; k++) {
          // might have moved due to reallocations
          SparseArray<T> &colj = colsB.getColumn (jindex);
          size_t kindex = rs.pMinus.keyAt (k);
          SparseArray<T> &colk = colsB.getColumn (kindex);
          // operate a linear combination on the columns of indices j and k
          // in order to get a new column having the pair.getFirst element equal
          // to zero
          T alpha = -colk.get (targetRow);
          T beta = colj.get (targetRow);
          T gcdt = std::gcd (alpha, beta);
          alpha /= gcdt;
          beta /= gcdt;

          if (DEBUG) {
            std::cout << "Computing " << beta << " * " << colk << " + " << alpha
                << " * " << colj << std::endl;
          }
          // single pivot forced by pure pos, or last pivot to be applied
          // modify the column colk in place
          if (purePos != -1 || j == je - 1) {
            auto changed = sumProdInto (beta, colk, alpha, colj);
            for (size_t ind = 0, inde = changed.size (); ind < inde; ind++) {
              size_t key = changed[ind].first;
              rowSigns.setValue (key, kindex, changed[ind].second);
            }
            normalize (colk);
            if (DEBUG) {
              std::cout << "Built (emplaced at " << kindex << ") " << colk
                  << std::endl;
              std::cout << "Changed reported " << changed << std::endl;
            }
          } else {
            auto newCol = sumProd (beta, colk, alpha, colj);
            normalize (newCol);
            if (newCol.size () > 0) {
              size_t newColIndex = colsB.getColumnCount ();
              for (size_t ind = 0, inde = newCol.size (); ind < inde; ind++) {
                rowSigns.setValue (newCol.keyAt (ind), newColIndex,
                                   newCol.valueAt (ind));
              }
              if (DEBUG) {
                std::cout << "Built (appended at " << newColIndex << ") "
                    << newCol << std::endl;
              }
              colsB.appendColumn (std::move (newCol));
            }
          }
          if (DEBUG) {
            std::cout << "Obtained rowsign " << rowSigns.get (targetRow)
                << " for row " << targetRow << std::endl;
          }
        }

      }

    }

    static void applyRowElimination (MatrixCol<T> &matC, MatrixCol<T> &matB,
    RowSigns<T> &rowSigns,
                                     std::pair<size_t, size_t> &counts,
                                     bool onlyPositive,
                                     const EliminationHeuristic &heur)
    {
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
      T cHk = matC.get (tRow, tCol);
      T bbeta = std::abs (cHk);
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
        // that means it exists a row that all components are positive respectivly negativ
        // [1.1.a] delete from the extended matrix all the columns of index j \in P+ \cup P-

        for (size_t i = 0, ie = toVisit.size (); i < ie; i++) {
          size_t j = toVisit.keyAt (i);
          matC.getColumn (j).clear ();
          matB.getColumn (j).clear ();
        }
        rowSigns.clearRow (tRow);
        return;
      }

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
        int tCol = isNeg ? rowData.pMinus.keyAt (0) : rowData.pPlus.keyAt (0);

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

    static PivotChoice findBestPivot (const MatrixCol<T> &matC,
                                      const RowSigns<T> &rowSigns,
                                      bool onlyPositive, size_t loopLimit)
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

    static ssize_t findBestFMERow (const MatrixCol<T> &matB,
                                   const RowSigns<T> &rowSigns,
                                   size_t loopLimit)
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
        int weight = pps * ppm - pps - ppm;

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

  public:
    static void clearColumn (int tCol, MatrixCol<T> &matC, MatrixCol<T> &matB,
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

  };

}

#endif /* INVARIANTCALCULATOR_H_ */
