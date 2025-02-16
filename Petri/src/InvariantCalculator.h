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

    class RowSign
    {

    public:
      // The row
      const size_t row;
      // P+ set
      SparseBoolArray pPlus;
      // P- set
      SparseBoolArray pMinus;

      /**
       * initially empty.
       *
       * @param row
       */
      RowSign (size_t row)
          : row (row), pPlus (), pMinus ()
      {
      }

      void setValue (size_t j, T val)
      {
        if (val == 0) {
          pMinus.clear (j);
          pPlus.clear (j);
        } else if (val < 0) {
          pMinus.set (j);
          pPlus.clear (j);
        } else {
          pMinus.clear (j);
          pPlus.set (j);
        }
      }

      void print (std::ostream &os) const
      {
        os << "RowSign [row=" << row << ", pPlus=" << pPlus << ", pMinus="
            << pMinus << "]";
        return;
      }

      friend std::ostream& operator<< (std::ostream &os, const RowSign &obj)
      {
        obj.print (os);
        return os;
      }

    };

    class RowSigns
    {
    public:
      // Nested class to hold sign information for a given row.
      /**
       * A class for holding the sets P+ = {j | c_hj &gt; 0} and P- = {j | c_hj less
       * 0} for a given row.
       */

      RowSigns (const MatrixCol<T> &matC)
      {
        size_t rowCount = matC.getRowCount ();
        // Preallocate one RowSign per row.
        rows.reserve (rowCount);
        for (size_t row = 0; row < rowCount; ++row) {
          rows.push_back (RowSign (row));
        }
        // For each column in matC, update the corresponding row.
        for (size_t icol = 0, cole = matC.getColumnCount (); icol < cole;
            ++icol) {
          const SparseArray<T> &col = matC.getColumn (icol);
          for (size_t i = 0, ie = col.size (); i < ie; ++i) {
            size_t row = col.keyAt (i);
            // Update row 'row' with the sign information from column icol.
            if (col.valueAt (i) < 0) {
              rows[row].pMinus.append (icol, true);
            } else {
              rows[row].pPlus.append (icol, true);
            }
          }
        }
      }
      // updateValue forces all modifications to go through RowSigns.
      // This way, any extra bookkeeping we need later is performed here.
      void setValue (size_t row, size_t col, T newVal)
      {
        // For now, simply delegate to the underlying RowSign.
        rows[row].setValue (col, newVal);
        // (Additional bookkeeping could be added here later.)
      }

      const RowSign& get (size_t row) const
      {
        return rows[row];
      }

      size_t size () const
      {
        return rows.size ();
      }

      ssize_t findSingleSignRow (size_t startIndex) const
      {
        size_t sz = rows.size ();
        // First, scan from startIndex to end.
        for (size_t i = startIndex; i < sz; ++i) {
          if (rows[i].pMinus.size () == 1 || rows[i].pPlus.size () == 1) return static_cast<ssize_t> (i);
        }
        // Then, scan from beginning up to startIndex.
        for (size_t i = 0; i < static_cast<size_t> (startIndex); ++i) {
          if (rows[i].pMinus.size () == 1 || rows[i].pPlus.size () == 1) return static_cast<ssize_t> (i);
        }
        return -1;
      }

    private:
      // The underlying container for row sign data.
      std::vector<RowSign> rows;

    };

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
        MatrixCol<T> mat, bool onlyPositive)
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

      MatrixCol<T> matB = phase1PIPE (matnorm.transpose ());

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
      while (colsB.getColumnCount () < 20000) {
        if (treated.size () > 0) {
          for (ssize_t i = treated.size () - 1; i >= 0; i--) {
            colsBsparse.insert (colsB.getColumn (treated.keyAt (i)));
            colsB.deleteColumn (treated.keyAt (i));
          }
          treated.clear ();
        }

        RowSigns rowSigns (colsB);
        SparseBoolArray negRows;

        int minRow = -1;
        int minRowWeight = -1;
        for (size_t row = 0, rowe = rowSigns.size (); row < rowe; row++) {
          RowSign pp = rowSigns.get (row);
          int pps = pp.pPlus.size ();
          int ppm = pp.pMinus.size ();
          int weight = pps + ppm;

          if (pps == 0) {
            for (int i = 0; i < ppm; i++) {
              negRows.set (pp.pMinus.keyAt (i));
            }
          }
          if (pps > 0 && ppm > 0) {
            if (pps == 1 || ppm == 1) {
              // can't grow the size
              minRow = row;
              break;
            }
            if (minRow == -1 || minRowWeight > weight) {
              int refinedweight = 0;
              for (size_t i = 0, ie = pp.pPlus.size (); i < ie; i++) {
                refinedweight += colsB.getColumn (pp.pPlus.keyAt (i)).size ();
              }
              for (size_t i = 0, ie = pp.pMinus.size (); i < ie; i++) {
                refinedweight += colsB.getColumn (pp.pMinus.keyAt (i)).size ();
              }
              if (minRow == -1 || minRowWeight > refinedweight) {
                minRow = row;
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
    static void removeNegativeValues (
        std::unordered_set<SparseArray<T>> &colsBsparse)
    {
      for (auto it = colsBsparse.begin (); it != colsBsparse.end ();) {
        const SparseArray<T> &a = *it;
        bool hasNegativeValue = false;
        for (size_t i = 0, ie = a.size (); i < ie; i++) {
          if (a.valueAt (i) < 0) {
            hasNegativeValue = true;
            break;
          }
        }
        if (hasNegativeValue) {
          it = colsBsparse.erase (it);
        } else {
          ++it;
        }
      }
    }

    /**
     * Efficiently removes trivial columns (i.e. empty columns in matC) from both
     * matC and matB. For each trivial column, its corresponding invariant from
     * matB is moved into trivialInv.
     *
     * @param matC       The incidence matrix.
     * @param matB       The transformation matrix.
     * @param trivialInv Vector accumulating trivial invariants.
     */
    static void cullConstantColumns (MatrixCol<T> &matC, MatrixCol<T> &matB,
                                     std::vector<SparseArray<T>> &trivialInv)
    {
      std::vector<size_t> trivialIndexes;
      size_t colCount = matC.getColumnCount ();

      // First pass: accumulate indexes of trivial columns.
      for (size_t col = 0; col < colCount; ++col) {
        if (matC.getColumn (col).size () == 0) {
          trivialIndexes.push_back (col);
        }
      }

      // Only rebuild if we found any trivial columns.
      if (!trivialIndexes.empty ()) {
        std::cout << "Culling trivial invariants: removed "
            << trivialIndexes.size () << " columns." << std::endl;

        // Prepare new vectors for non-trivial columns.
        std::vector<SparseArray<T>> newColsC;
        std::vector<SparseArray<T>> newColsB;
        newColsC.reserve (colCount - trivialIndexes.size ());
        newColsB.reserve (colCount - trivialIndexes.size ());

        size_t trivialPos = 0; // Pointer into trivialIndexes (sorted in increasing order)
        for (size_t col = 0; col < colCount; ++col) {
          if (trivialPos < trivialIndexes.size ()
              && col == trivialIndexes[trivialPos]) {
            // Found a trivial column: move its corresponding invariant from matB.
            trivialInv.push_back (std::move (matB.getColumn (col)));
            ++trivialPos;
            // Skip this column from matC.
          } else {
            // Non-trivial column: move from both matrices.
            newColsC.push_back (std::move (matC.getColumn (col)));
            newColsB.push_back (std::move (matB.getColumn (col)));
          }
        }
        // Replace the old columns with the new (culled) ones.
        matC.setColumns (std::move (newColsC));
        matB.setColumns (std::move (newColsB));
      }
    }

    /**
     * Efficiently removes duplicate columns in matC.
     *
     * For each duplicate column k (i.e. a column equal to a previously seen column r),
     * we compute the difference diff = matB(k) - matB(r) using SparseArray::sumProd,
     * add it to trivialInv, and record k as a duplicate.
     *
     * After processing all columns, if any duplicates were found, we rebuild matC and matB
     * (skipping the duplicate columns) in one sweep.
     *
     * @tparam T            Integral type.
     * @param matC          The incidence matrix.
     * @param matB          The transformation matrix.
     * @param trivialInv    Vector accumulating trivial invariants.
     */
    static void cullDuplicateColumns (MatrixCol<T> &matC, MatrixCol<T> &matB,
                                      std::vector<SparseArray<T>> &trivialInv)
    {
      size_t colCount = matC.getColumnCount ();
      std::vector<size_t> duplicateIndexes;  // collect duplicate column indices
      std::unordered_map<SparseArray<T>, size_t> repMap;

      // First pass: detect duplicates.
      for (size_t col = 0; col < colCount; ++col) {
        const SparseArray<T> &columnC = matC.getColumn (col);
        auto it = repMap.find (columnC);
        if (it == repMap.end ()) {
          repMap.insert (
            { columnC, col });
        } else {
          // Duplicate column found: record its index.
          duplicateIndexes.push_back (col);
          // Compute diff = matB(col) - matB(rep) using sumProd.
          SparseArray<T> diff = SparseArray<T>::sumProd (
              1, matB.getColumn (col), -1, matB.getColumn (it->second));
          trivialInv.push_back (std::move (diff));
        }
      }

      // Only rebuild if we found duplicates.
      if (!duplicateIndexes.empty ()) {
        std::cout << "Culling duplicate invariants: removed "
            << duplicateIndexes.size () << " columns." << std::endl;

        // duplicateIndexes are produced in increasing order since col increases.
        std::vector<SparseArray<T>> newColsC;
        std::vector<SparseArray<T>> newColsB;
        newColsC.reserve (colCount - duplicateIndexes.size ());
        newColsB.reserve (colCount - duplicateIndexes.size ());

        size_t dupPos = 0; // pointer into duplicateIndexes
        for (size_t col = 0; col < colCount; ++col) {
          if (dupPos < duplicateIndexes.size ()
              && duplicateIndexes[dupPos] == col) {
            ++dupPos; // skip this duplicate column
          } else {
            newColsC.push_back (std::move (matC.getColumn (col)));
            newColsB.push_back (std::move (matB.getColumn (col)));
          }
        }

        // Replace the old columns with the new (culled) ones.
        matC.setColumns (std::move (newColsC));
        matB.setColumns (std::move (newColsB));
      }
    }

    static MatrixCol<T> phase1PIPE (MatrixCol<T> matC)
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
          << matC.getColumnCount () << " cols" << std::endl;
      RowSigns rowSigns (matC);
      int startIndex = 0;
      while (!matC.isZero ()) {
        startIndex = applyRowElimination (matC, matB, rowSigns, startIndex);
        if (DEBUG) {
          std::cout << "Mat max : " << matC.maxVal () << std::endl;
          std::cout << "B max : " << matB.maxVal () << std::endl;
        }
      }

      // Re-add the trivial invariants back into matB.
      for (const auto &inv : trivialInv) {
        matB.appendColumn (inv);
      }

      return matB;
    }

    static int applyRowElimination (MatrixCol<T> &matC, MatrixCol<T> &matB, RowSigns &rowSigns,
                       int startIndex)
    {
      // Find the candidate row with a single sign entry.
      ssize_t candidateRow = rowSigns.findSingleSignRow (startIndex);
      if (candidateRow != -1) {
        // Use candidateRow (cast to int if necessary) in test1b1.
        applySingleSignRowElimination (matC, matB, rowSigns, static_cast<size_t> (candidateRow));
        startIndex = candidateRow;
      } else {
        applyGeneralRowElimination (matC, matB, rowSigns);
      }
      return startIndex;
    }

  private:

    static void applyGeneralRowElimination (MatrixCol<T> &matC, MatrixCol<T> &matB,
                         RowSigns &rowSigns)
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

      T cHk = matC.get (tRow, tCol);
      T bbeta = std::abs (cHk);

      if (DEBUG) {
        std::cout << "Rule 1b2 : " << tCol << std::endl;
      }
      // for all cols j with j != tCol and c[tRow][j] != 0
      const RowSign &rowppm = rowSigns.get (tRow);
      if (DEBUG) {
        std::cout << "tCol : " << tCol << " tRow " << tRow << std::endl;
        std::cout << "rowppm : " << rowppm << std::endl;
      }
      SparseBoolArray toVisit = SparseBoolArray::unionOperation (rowppm.pMinus,
                                                                 rowppm.pPlus);

      for (size_t i = 0; i < toVisit.size (); i++) {
        size_t j = toVisit.keyAt (i);
        SparseArray<T> &colj = matC.getColumn (j);

        if (j == tCol) {
          continue;
        }

        T cHj = colj.get (tRow);
        if (cHj != 0) {
          // substitute to the column of index j the linear combination
          // of the columns of indices tCol and j with coefficients
          // alpha and beta defined as follows:
          T alpha =
              ((signum (cHj) * signum (cHk)) < 0) ?
                  std::abs (cHj) : -std::abs (cHj);
          if (alpha == 0 && bbeta == 1) {
            continue;
          }
          T gcdt = std::gcd (alpha, bbeta);
          alpha /= gcdt;
          T beta = bbeta / gcdt;

          if (DEBUG) {
            std::cout << "rule 1. b 2 : " << alpha << "*" << tCol << " + "
                << beta << " * " << j << std::endl;
            std::cout << "tCol : " << matC.getColumn (tCol) << std::endl;
            std::cout << "colj : " << colj << std::endl;
          }

          SparseBoolArray changed = sumProdInto (beta, colj, alpha,
                                                 matC.getColumn (tCol));
          for (size_t ind = 0, inde = changed.size (); ind < inde; ind++) {
            rowSigns.setValue (changed.keyAt (ind), j,
                            colj.get (changed.keyAt (ind)));
          }
          SparseArray<T> &coljb = matB.getColumn (j);
          if (DEBUG) {
            std::cout << "colj(after) : " << colj << std::endl;
            std::cout << "B[colj] : " << coljb << std::endl;
          }

          sumProdInto (beta, coljb, alpha, matB.getColumn (tCol));

          T gcdm = gcd (matC.getColumn (j));
          if (gcdm != 1) {
            T gcdb = gcd (matB.getColumn (j));
            if (gcdb != 1) {
              T ggcd = std::gcd (gcdm, gcdb);
              if (ggcd != 1) {
                matC.getColumn (j).scalarDiv (ggcd);
                matB.getColumn (j).scalarDiv (ggcd);
                //std::cout << "1b1 reduce by " << ggcd << std::endl;
              }
            }
          }

          if (DEBUG) {
            std::cout << "B[tCol] : " << matB.getColumn (tCol) << std::endl;
            std::cout << " after B[colj] : " << coljb << std::endl;
          }
        }
      }
      clearColumn (tCol, matC, matB, rowSigns);
    }

  public:
    static void clearColumn (int tCol, MatrixCol<T> &matC, MatrixCol<T> &matB,
                             RowSigns &rowSigns)
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
                         RowSigns &rowSigns, size_t candidateRow)
    {
      if (DEBUG) {
        std::cout << "Rule 1b.1 : " << candidateRow << std::endl;
      }
      // Get the candidate row data freshly.
      const RowSign &rowData = rowSigns.get (candidateRow);
      // In our construction, exactly one of pPlus or pMinus must have size 1.
      assert(rowData.pPlus.size() == 1 || rowData.pMinus.size() == 1);

      // Determine which set is unique.
      // If pPlus is unique then tCol is the unique key and the complementary set is pMinus;
      // otherwise, tCol is from pMinus and the complementary set is pPlus.
      bool isPos = (rowData.pPlus.size () == 1);
      int tCol = isPos ? rowData.pPlus.keyAt (0) : rowData.pMinus.keyAt (0);

      // Loop while the complementary set (re-read fresh each iteration) is non-empty.
      while (true) {
        // Re-read the candidate row to get the latest state.
        const RowSign &currentRow = rowSigns.get (candidateRow);
        const SparseBoolArray &currentComplement =
            isPos ? currentRow.pMinus : currentRow.pPlus;
        if (currentComplement.size () == 0) {
          break;
        }
        int j = currentComplement.keyAt (0);

        // Retrieve the coefficients from the candidate row.
        T chk = std::abs (matC.get (candidateRow, tCol));
        T chj = std::abs (matC.get (candidateRow, j));
        T gcdt = std::gcd (chk, chj);
        chk /= gcdt;
        chj /= gcdt;

        // Update matC: combine columns j and tCol.
        SparseBoolArray changed = sumProdInto (chk, matC.getColumn (j), chj,
                                               matC.getColumn (tCol));
        // For each change, update the row-sign bookkeeping.
        // (We re-read each row via rowSigns.setValue to avoid caching any references.)
        for (size_t ind = 0, inde = changed.size (); ind < inde; ++ind) {
          size_t key = changed.keyAt (ind);
          rowSigns.setValue (key, j, matC.getColumn (j).get (key));
        }

        // Update matB with the same linear combination.
        SparseArray<T> &coljb = matB.getColumn (j);
        sumProdInto (chk, coljb, chj, matB.getColumn (tCol));

        // Optionally perform a GCD reduction on column j.
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
        // Loop condition re-evaluated by re-reading rowSigns.get(candidateRow)
      }
      // Finally, clear the candidate column from both matrices.
      clearColumn (tCol, matC, matB, rowSigns);
    }


  };

}

#endif /* INVARIANTCALCULATOR_H_ */
