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

namespace petri {

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
        MatrixCol<T> mat, bool onlyPositive,  const EliminationHeuristic &heur=EliminationHeuristic())
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

      MatrixCol<T> matB = phase1PIPE (matnorm.transpose (), heur);

//		const MatrixCol<T> matB = phase1PIPE(new MatrixCol<T>(mat));
      // We want to work with columns in this part of the algorithm
      // We add and remove columns all day => we want to switch to a column based
      // representation
      // order of rows is really irrelevant + columns which are identical up to
      // scaling factor are useless
      // let's use a set of columns.
      std::unordered_set<SparseArray<T>> colsBsparse (
          2 * matB.getColumnCount ());
      for (size_t i = 0; i < matB.getColumnCount (); i++) {
        SparseArray<T> &col = matB.getColumn (i);
        if (col.size () != 0) {
          normalizeWithSign (col);
          colsBsparse.insert (col);
        }
      }

      if (!onlyPositive) {
        return colsBsparse;
      }

      MatrixCol<T> colsB (tmat.getRowCount (), 0);
      for (const SparseArray<T> &cb : colsBsparse) {
        colsB.appendColumn (cb);
      }

      // phase 2
      std::cout << "// Phase 2 : computing semi flows from basis of "
          << colsB.getColumnCount () << " invariants " << std::endl;

      //int iter = 0;
      SparseBoolArray treated;
      colsBsparse = std::unordered_set<SparseArray<T>> ();
      while (true || colsB.getColumnCount () < 20000) {
        if (treated.size () > 0) {
          for (ssize_t i = treated.size () - 1; i >= 0; i--) {
            colsBsparse.insert (colsB.getColumn (treated.keyAt (i)));
            colsB.deleteColumn (treated.keyAt (i));
          }
          treated.clear ();
        }

        RowSigns<T> rowSigns (colsB);
        SparseBoolArray negRows;

        int minRow = -1;
        int minRowWeight = -1;
        for (const auto &rs : rowSigns) {
          int pps = rs.pPlus.size();
          int ppm = rs.pMinus.size();
          int weight = pps + ppm;

          if (pps == 0) {
            for (int i = 0; i < ppm; i++) {
              negRows.set(rs.pMinus.keyAt(i));
            }
          }
          if (pps > 0 && ppm > 0) {
            if (pps == 1 || ppm == 1) {
              // can't grow the size; use the stored row index
              minRow = rs.row;
              break;
            }
            if (minRow == -1 || minRowWeight > weight) {
              int refinedweight = 0;
              for (size_t i = 0, ie = rs.pPlus.size(); i < ie; i++) {
                refinedweight += colsB.getColumn(rs.pPlus.keyAt(i)).size();
              }
              for (size_t i = 0, ie = rs.pMinus.size(); i < ie; i++) {
                refinedweight += colsB.getColumn(rs.pMinus.keyAt(i)).size();
              }
              if (minRow == -1 || minRowWeight > refinedweight) {
                minRow = rs.row;
                minRowWeight = refinedweight;
              }
            }
          }
        }

        if (negRows.size () > 0) {
          // cleanup
          for (ssize_t j = negRows.size () - 1; j >= 0; j--) {
            colsB.deleteColumn (negRows.keyAt (j));
          }
          continue;
        }
        // check for a pure positive column
        int purePos = -1;

        for (size_t i = 0, ie = colsB.getColumnCount (); i < ie; i++) {
          if (treated.get (i)) {
            continue;
          }
          SparseArray<T> &col = colsB.getColumn (i);
          bool hasNeg = false;
          for (size_t j = 0, je = col.size (); j < je; j++) {
            if (col.valueAt (j) < 0) {
              hasNeg = true;
              break;
            }
          }
          if (!hasNeg) {
            // check intersection
            bool needed = false;
            for (size_t j = 0, je = col.size (); j < je; j++) {
              int row = col.keyAt (j);
              RowSign ppm = rowSigns.get (row);
              if (ppm.pMinus.size () > 0) {
                needed = true;
                purePos = i;
                minRow = row;
                break;
              }
            }
            if (!needed) {
              treated.set (i);
            } else {
              break;
            }
          }
        }

        int targetRow = minRow;
        if (targetRow == -1) {
          // no more negative rows to treat
          break;
        }
        RowSign ppm = rowSigns.get (targetRow);
        if (ppm.pPlus.size () > 0) {
          for (size_t j = 0, je = ppm.pPlus.size (); j < je; j++) {
            auto jindex = ppm.pPlus.keyAt (j);
            if (purePos != -1) {
              jindex = purePos;
              j = je;
            }
            for (size_t k = 0, ke = ppm.pMinus.size (); k < ke; k++) {
              // might have moved due to reallocations
              SparseArray<T> &colj = colsB.getColumn (jindex);
              SparseArray<T> &colk = colsB.getColumn (ppm.pMinus.keyAt (k));
              // operate a linear combination on the columns of indices j and k
              // in order to get a new column having the pair.getFirst element equal
              // to zero
              int a = -colk.get (targetRow);
              int b = colj.get (targetRow);
              SparseArray<T> column = SparseArray<T>::sumProd (a, colj, b,
                                                               colk);
              // add normalization step : we don't need scalar scaling of each other
              normalize (column);
              // append column to matrix B
              // tests existence
              if (column.size () > 0) {
                colsB.appendColumn (column);
              }
            }
          }
          // Delete from B all the columns of index k \in P-
          // cleanup
          for (ssize_t j = ppm.pMinus.size () - 1; j >= 0; j--) {
            colsB.deleteColumn (ppm.pMinus.keyAt (j));
            treated.deleteAndShift (ppm.pMinus.keyAt (j));
          }
        }
        // std::cout << "Phase 2 iter " << iter++ << " rows : " <<
        // colsB.getRowCount() << " cols " << colsB.getColumnCount() << " treated " <<
        // colsBsparse.size() << std::endl;
        // std::cout << colsB << std::endl;
      }
      // std::cout << "Found "<< colsB.getColumnCount() << " invariants."<< std::endl;

      for (SparseArray<T> l : colsB.getColumns ()) {
        if (l.size () > 0) {
          colsBsparse.insert (l);
        }
      }
      // std::cout << "Found "<< colsBsparse.size() << " different invariants."<< std:endl;
      removeNegativeValues (colsBsparse);
      std::cout << "Found " << colsBsparse.size () << " positive invariants."
          << std::endl;
      return colsBsparse;
    }

  private:


    static MatrixCol<T> phase1PIPE (MatrixCol<T> matC, const EliminationHeuristic &heur)
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
      size_t nbArcs =0;
      for (const auto &col : matC.getColumns ()) {
        nbArcs += col.size ();
      }
      std::cout << "// Phase 1: matrix " << matC.getRowCount () << " rows "
          << matC.getColumnCount () << " cols " << nbArcs << " entries" << std::endl;
      RowSigns rowSigns (matC,heur.useSingleSignRow());

      std::pair<size_t,size_t> counts (0,0);
      int startIndex = 0;
      while (!matC.isZero ()) {
        startIndex = applyRowElimination (matC, matB, rowSigns, startIndex, counts, heur);
        if (DEBUG) {
          std::cout << "Mat max : " << matC.maxVal () << std::endl;
          std::cout << "B max : " << matB.maxVal () << std::endl;
        }
      }
      std::cout << "Finished phase 1 with " << counts.first << " SingleSign rule and " << counts.second << " generalized " << std::endl;
      // Re-add the trivial invariants back into matB.
      for (const auto &inv : trivialInv) {
        matB.appendColumn (inv);
      }

      return matB;
    }

    static int applyRowElimination(MatrixCol<T>& matC,
                                   MatrixCol<T>& matB,
                                   RowSigns<T>& rowSigns,
                                   int startIndex,
                                   std::pair<size_t,size_t> &counts, const EliminationHeuristic &heur)
    {
        // 1) Check Single-Sign rows first:
      if (heur.useSingleSignRow()) {
        ssize_t candidateRow = rowSigns.findSingleSignRow(startIndex, heur.getLoopLimit());
        if (candidateRow != -1) {
            // Single-sign path => pick pivot col from pMinus or pPlus
            // possibly compare matC column sizes if both are size 1
            applySingleSignRowElimination(matC, matB, rowSigns, static_cast<size_t>(candidateRow));
            counts.first++;
            return candidateRow;
        }
      }
        // 2) General pivot choice:
        //    (a) pick a column with minimal size in matC
        //    (b) among the rows in that column, pick tRow with the smallest "cost"
        //        (for instance, the smallest absolute cell value, or minimal expansions)
        auto pivot = findBestPivot(matC, rowSigns, heur.getLoopLimit());
        // pivot is a struct { size_t row; size_t col; }
        eliminateRowWithPivot(pivot.row, pivot.col, matC, matB, rowSigns);
        counts.second++;
        return startIndex;  // or possibly pivot.row, depending on your logic
    }


    static int applyRowElimination2 (MatrixCol<T> &matC, MatrixCol<T> &matB, RowSigns<T> &rowSigns,
                       int startIndex, std::pair<size_t,size_t> & counts)
    {
      // Find the candidate row with a single sign entry.
      ssize_t candidateRow = rowSigns.findSingleSignRow (startIndex);
      if (candidateRow != -1) {
        // Use candidateRow (cast to int if necessary) in test1b1.
        applySingleSignRowElimination (matC, matB, rowSigns, static_cast<size_t> (candidateRow));
        startIndex = candidateRow;
        counts.first++;
      } else {
        applyGeneralRowElimination (matC, matB, rowSigns);
        counts.second++;
      }
      return startIndex;
    }

  private:

    static void eliminateRowWithPivot(size_t tRow, int tCol,
                                      MatrixCol<T>& matC,
                                      MatrixCol<T>& matB,
                                      RowSigns<T>& rowSigns)
    {
        T cHk = matC.get(tRow, tCol);
        T bbeta = std::abs(cHk);
        const auto & rowsign = rowSigns.get(tRow);
        SparseBoolArray toVisit = SparseBoolArray::unionOperation (rowsign.pMinus,
                                                                         rowsign.pPlus);

        if (DEBUG) {
          std::cout << "tCol : " << tCol << " tRow " << tRow << std::endl;
          std::cout << "rowsign : " << rowsign << std::endl;
        }

        for (size_t i = 0; i < toVisit.size(); i++) {
            size_t j = toVisit.keyAt(i);
            if (j == (size_t)tCol) {
                continue;
            }

            SparseArray<T>& colj = matC.getColumn(j);
            T cHj = colj.get(tRow);
            if (cHj != 0) {
                // Compute alpha and beta
                T alpha = ((signum(cHj) * signum(cHk)) < 0)
                          ? std::abs(cHj)
                          : -std::abs(cHj);
                if (alpha == 0 && bbeta == 1) {
                    // No operation needed
                    continue;
                }
                T gcdt = std::gcd(alpha, bbeta);
                alpha /= gcdt;
                T beta = bbeta / gcdt;

                if (DEBUG) {
                    std::cout << "Eliminating with pivot col="
                              << tCol << " row=" << tRow
                              << " combining j=" << j << std::endl;
                }

                // Update matC
                SparseArray<T> changed = sumProdInto(beta,
                                                      colj,
                                                      alpha,
                                                      matC.getColumn(tCol));
                // Update the row-sign bookkeeping
                for (size_t ind = 0, inde = changed.size(); ind < inde; ind++) {
                    size_t key = changed.keyAt(ind);
                    rowSigns.setValue(key, j, changed.valueAt(ind));
                }

                // Update matB
                SparseArray<T>& coljb = matB.getColumn(j);
                sumProdIntoNoChange(beta, coljb, alpha, matB.getColumn(tCol));

                // Optionally perform GCD reduction on column j
                T gcdm = gcd(matC.getColumn(j));
                if (gcdm != 1) {
                    T gcdb = gcd(matB.getColumn(j));
                    if (gcdb != 1) {
                        T ggcd = std::gcd(gcdm, gcdb);
                        if (ggcd != 1) {
                            matC.getColumn(j).scalarDiv(ggcd);
                            matB.getColumn(j).scalarDiv(ggcd);
                        }
                    }
                }
            }
        }
        // Finally clear the pivot column tCol
        clearColumn(tCol, matC, matB, rowSigns);
    }

    struct PivotChoice {
        size_t row;
        size_t col;
    };

    static PivotChoice findBestPivot(const MatrixCol<T>& matC,
                                     const RowSigns<T>& rowSigns, size_t loopLimit)
    {
        // We'll assume there's at least one non-empty column; otherwise we can't pivot.
        // If you want to handle an all-empty matrix, you can do so by checking further.

        size_t bestColSize = std::numeric_limits<size_t>::max();
        std::vector<size_t> candidateCols;
        candidateCols.reserve(matC.getColumnCount());

        // 1) Find the columns that have the minimal nonzero size in one pass.
        for (size_t c = 0; c < matC.getColumnCount(); c++) {
            size_t s = matC.getColumn(c).size();
            if (s == 0) {
                continue; // skip empty columns
            }
            if (s < bestColSize) {
                bestColSize = s;
                candidateCols.clear();
                candidateCols.push_back(c);
            }
            else if (s == bestColSize) {
                candidateCols.push_back(c);
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
        size_t bestRow    = 0;
        size_t bestCol    = candidateCols[0];       // default
        size_t bestRowSz  = std::numeric_limits<size_t>::max();
        T      bestAbsVal = std::numeric_limits<T>::max();
        bool   foundPivot = false;

        // For each candidate column, iterate over its nonzero entries
        for (size_t c: candidateCols) {
            const SparseArray<T>& colData = matC.getColumn(c);
            for (size_t i = 0, ie = colData.size(); i < ie; i++) {
                size_t r = colData.keyAt(i);
                T val    = colData.valueAt(i);
                // We assume val != 0 in a properly stored SparseArray.

                // Row size in rowSigns
                const auto& rs = rowSigns.get(r);
                size_t rowSz   = rs.pPlus.size() + rs.pMinus.size();
                T absVal       = (val < 0) ? -val : val;

                // Compare to see if it's better:
                //   - pick the pivot that has the smaller rowSz
                //   - tie-break on smaller absVal
                if (!foundPivot ||
                    (rowSz < bestRowSz) ||
                    (rowSz == bestRowSz && absVal < bestAbsVal))
                {
                    foundPivot  = true;
                    bestRow     = r;
                    bestCol     = c;
                    bestRowSz   = rowSz;
                    bestAbsVal  = absVal;

                    // Early‐exit if (rowSz <= 2) and (absVal == 1).
                    if (rowSz <= 2 && absVal == 1) {
                        // Immediately return this pivot.
                        return PivotChoice { bestRow, bestCol };
                    }
                }
            }
        }

        // If we found no pivot at all, that means all columns might have size 0
        // or we never had nonzero entries. But in practice we do handle above.
        // So we finalize the best we found:
        return foundPivot
            ? PivotChoice { bestRow, bestCol }
            : PivotChoice { 0, 0 };  // or throw/assert if truly impossible
    }



    static void applyGeneralRowElimination (MatrixCol<T> &matC, MatrixCol<T> &matB,
                         RowSigns<T> &rowSigns)
    {
      // [1.1.b.1] let tRow be the index of a non-zero row of C.
      // let tCol be the index of a column such that c[trow][tcol] != 0.

      ssize_t candidate = -1;
      size_t szcand = std::numeric_limits<T>::max ();
      size_t totalcand = std::numeric_limits<T>::max ();
      for (size_t col = 0; col < matC.getColumnCount (); col++) {
        size_t size = matC.getColumn (col).size ();
        if (size == 0) {
          continue;
        } else if (size <= szcand) {
          size_t total = sumAbsValues (matC.getColumn (col));
          if (size < szcand || (size == szcand && total <= totalcand)) {
            candidate = col;
            szcand = size;
            totalcand = total;
          }
        }
      }
      size_t tRow = matC.getColumn (candidate).keyAt (0);
      size_t tCol = candidate;


      eliminateRowWithPivot(tRow, tCol, matC, matB, rowSigns);
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

  private:

    static void applySingleSignRowElimination (MatrixCol<T> &matC, MatrixCol<T> &matB,
                         RowSigns<T> &rowSigns, size_t tRow)
    {
      if (DEBUG) {
        std::cout << "Rule 1b.1 : " << tRow << std::endl;
      }
      // Get the candidate row data freshly.
      const auto &rowData = rowSigns.get (tRow);
      // In our construction, exactly one of pPlus or pMinus must have size 1.
      assert(rowData.pPlus.size() == 1 || rowData.pMinus.size() == 1);

      // Determine which set is unique.
      // If pPlus is unique then tCol is the unique key and the complementary set is pMinus;
      // otherwise, tCol is from pMinus and the complementary set is pPlus.
      bool isNeg = (rowData.pMinus.size() == 1);
      int tCol = isNeg ? rowData.pMinus.keyAt(0) : rowData.pPlus.keyAt(0);

      if (rowData.pPlus.size () == 1 && rowData.pMinus.size () == 1) {
        // std::cout << "Examine col j=" << j << " size=" << matC.getColumn (j).size () <<" with col=" << tCol << " size " << matC.getColumn (tCol).size () << std::endl;
        // we can actually choose which one to get rid of
        if (matC.getColumn (rowData.pMinus.keyAt(0)).size () > matC.getColumn (rowData.pPlus.keyAt(0)).size ()) {
          tCol = rowData.pPlus.keyAt(0);
          isNeg = !isNeg;
          // std::cout << "SWAPPED" << std::endl;
        }
      }

      eliminateRowWithPivot(tRow, tCol, matC, matB, rowSigns);
    }


  };

}

#endif /* INVARIANTCALCULATOR_H_ */
