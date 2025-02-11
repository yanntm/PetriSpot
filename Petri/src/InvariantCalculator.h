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
     * A class for holding the sets P+ = {j | c_hj &gt; 0} and P- = {j | c_hj less
     * 0} for a given row.
     */
    class PpPm
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
      PpPm (size_t row)
          : row (row), pPlus (), pMinus ()
      {
      }

      void setValue (size_t j, int val)
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
        os << "PpPm [row=" << row << ", pPlus=" << pPlus << ", pMinus="
            << pMinus << "]";
        return;
      }

      friend std::ostream& operator<< (std::ostream &os, const PpPm &obj)
      {
        obj.print (os);
        return os;
      }

    };

    /**
     * Calculates for a given matrix the P+ = {j | c_hj &gt; 0} and P- = {j | c_hj
     * less 0} sets for each row.
     *
     * @param matC - the matrix from which the sets should be calculated.
     * @return The result of the calculation
     */
    static std::vector<PpPm> calcPpPm (const MatrixCol<T> &matC)
    {
      std::vector<PpPm> result;
      for (size_t row = 0; row < matC.getRowCount (); row++) {
        result.push_back (PpPm (row));
      }
      for (size_t icol = 0, cole = matC.getColumnCount (); icol < cole;
          icol++) {
        const SparseArray<T> &col = matC.getColumn (icol);
        for (size_t i = 0, ie = col.size (); i < ie; i++) {
          PpPm &toedit = result[col.keyAt (i)];
          if (col.valueAt (i) < 0) {
            toedit.pMinus.append (icol, true);
          } else {
            toedit.pPlus.append (icol, true);
          }
        }
      }
      return result;
    }

    class Check11bResult
    {

    public:
      // The first columnindex where c_hj < 0 respectivly c_hj > 0
      const int col;
      // The whole row
      const int row;
      // The set P+ respectivly P-
      SparseBoolArray *p;

      /**
       * Constructor to save the data.
       *
       * @paramk - the first column index where a component is less respectively
       *         greater than zero.
       * @param h     - the whole row.
       * @param pPlus - the set P+ respectively P-.
       */
      Check11bResult (int k, int row, SparseBoolArray *pPlus)
          : col (k), row (row), p (pPlus)
      {
      }
    };

    /**
     * Checks if there exists a row in the given matrix such that |P+| == 1 or |P-|
     * == 1 and returns the row or null if such a row do not exists.
     *
     * @param pppms      the list of all rows with P+ and P- sets.
     * @param startIndex
     * @return the row which satisfy |P+| == 1 or |P-| == 1 or null if not existent.
     */
    static Check11bResult check11b (std::vector<PpPm> &pppms, int startIndex)
    {

      for (auto it = pppms.begin () + startIndex; it != pppms.end (); ++it) {
        Check11bResult res = check11bPppm (*it);
        if (res.col != -1) {
          return res;
        }
      }
      for (auto it = pppms.begin (), _end = it + startIndex; it != _end; ++it) {
        Check11bResult res = check11bPppm (*it);
        if (res.col != -1) {
          return res;
        }
      }
      return Check11bResult (-1, -1, nullptr);
    }

    static Check11bResult check11bPppm (PpPm &pppm)
    {
      if (pppm.pMinus.size () == 1) {
        return Check11bResult (pppm.pMinus.keyAt (0), pppm.row, &pppm.pPlus);
      } else if (pppm.pPlus.size () == 1) {
        return Check11bResult (pppm.pPlus.keyAt (0), pppm.row, &pppm.pMinus);
      }
      return Check11bResult (-1, -1, nullptr);
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

        std::vector<PpPm> pppms = calcPpPm (colsB);
        SparseBoolArray negRows;

        int minRow = -1;
        int minRowWeight = -1;
        for (size_t row = 0, rowe = pppms.size (); row < rowe; row++) {
          PpPm pp = pppms[row];
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
              PpPm ppm = pppms[row];
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
        PpPm ppm = pppms[targetRow];
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

    static std::vector<int> normalize (SparseArray<T> &col, size_t size)
    {
      std::vector<int> list (size);
      bool allneg = true;
      for (size_t i = 0; i < col.size (); i++) {
        if (col.valueAt (i) > 0) {
          allneg = false;
          break;
        }
      }
      if (allneg) {
        for (size_t i = 0; i < col.size (); i++) {
          col.setValueAt (i, -col.valueAt (i));
        }
      }

      for (size_t i = 0; i < size; i++) {
        list.push_back (col.get (i, 0));
      }
      normalize (list);
      return list;
    }

    static MatrixCol<T> phase1PIPE (MatrixCol<T> matC)
    {
      // incidence matrix
      MatrixCol<T> matB = MatrixCol<T>::identity (matC.getColumnCount (),
                                                  matC.getColumnCount ());

      std::cout << "// Phase 1: matrix " << matC.getRowCount () << " rows "
          << matC.getColumnCount () << " cols" << std::endl;
      std::vector<PpPm> pppms = calcPpPm (matC);
      int startIndex = 0;
      while (!matC.isZero ()) {
        startIndex = test1b (matC, matB, pppms, startIndex);
        if (DEBUG) {
                std::cout << "Mat max : " << matC.maxVal() << std::endl;
                std::cout << "B max : " << matB.maxVal() << std::endl;
              }
      }
      return matB;
    }

    static int test1b (MatrixCol<T> &matC, MatrixCol<T> &matB,
                       std::vector<PpPm> &pppms, int startIndex)
    {
      // [1.1.b] if there exists a row h in C such that |P+| == 1 or |P-| == 1
      Check11bResult chkResult = check11b (pppms, startIndex);
      if (chkResult.col != -1) {
        test1b1 (matC, matB, pppms, chkResult);
        startIndex = chkResult.row;
      } else {
        test1b2 (matC, matB, pppms);
      }
      return startIndex;
    }

  public:
    static SparseBoolArray sumProdInto (int alpha, SparseArray<T> &ta, int beta,
                                        const SparseArray<T> &tb)
    {
      SparseBoolArray changed;
      SparseArray<T> flow (std::max (ta.size (), tb.size ()));

      size_t i = 0;
      size_t j = 0;
      while (i < ta.size () || j < tb.size ()) {
        unsigned int ki =
            i == ta.size () ?
                std::numeric_limits<unsigned int>::max () : ta.keyAt (i);
        unsigned int kj =
            j == tb.size () ?
                std::numeric_limits<unsigned int>::max () : tb.keyAt (j);
        if (ki == kj) {
          T val = petri::addExact (petri::multiplyExact (alpha, ta.valueAt (i)),
                                   petri::multiplyExact (beta, tb.valueAt (j)));
          if (val != 0) {
            flow.append (ki, val);
          }
          if (val != ta.valueAt (i)) {
            changed.set (ki);
          }
          i++;
          j++;
        } else if (ki < kj) {
          T val = petri::multiplyExact (alpha, ta.valueAt (i));
          if (val != 0) {
            flow.append (ki, val);
          }
          if (val != ta.valueAt (i)) {
            changed.set (ki);
          }
          i++;
        } else if (kj < ki) {
          T val = petri::multiplyExact (beta, tb.valueAt (j));
          if (val != 0) {
            flow.append (kj, val);
          }
          if (val != 0) {
            changed.set (kj);
          }
          j++;
        }
      }
      ta = std::move(flow);
      return changed;
    }

  private:
    static int signum (T value)
    {
      return (value > 0) - (value < 0);
    }

    static void test1b2 (MatrixCol<T> &matC, MatrixCol<T> &matB,
                         std::vector<PpPm> &pppms)
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
      PpPm & rowppm = pppms[tRow];
      if (DEBUG) {
        std::cout << "tCol : " << tCol <<  " tRow " << tRow <<  std::endl;
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
            std::cout << "rule 1. b 2 : " << alpha << "*" << tCol << " + " << beta << " * " << j << std::endl;
            std::cout << "tCol : " << matC.getColumn (tCol) << std::endl;
            std::cout << "colj : " << colj << std::endl;
          }

          SparseBoolArray changed = sumProdInto (beta, colj, alpha,
                                                 matC.getColumn (tCol));
          for (size_t ind = 0, inde = changed.size (); ind < inde; ind++) {
            pppms[changed.keyAt (ind)].setValue (
                j, colj.get (changed.keyAt (ind)));
          }
          SparseArray<T> &coljb = matB.getColumn (j);
          if (DEBUG) {
            std::cout << "colj(after) : " << colj << std::endl;
            std::cout << "B[colj] : " << coljb << std::endl;
          }

          sumProdInto (beta, coljb, alpha, matB.getColumn (tCol));

          T gcdm = gcd(matC.getColumn (j));
          if (gcdm != 1) {
            T gcdb = gcd(matB.getColumn (j));
            if (gcdb != 1) {
              T ggcd = std::gcd(gcdm,gcdb);
              if (ggcd != 1) {
                matC.getColumn(j).scalarDiv(ggcd);
                matB.getColumn (j).scalarDiv(ggcd);
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
      clearColumn (tCol, matC, matB, pppms);
    }

  public:
    static void clearColumn (int tCol, MatrixCol<T> &matC, MatrixCol<T> &matB,
                             std::vector<PpPm> &pppms)
    {
      // delete from the extended matrix the column of index k
      SparseArray<T> &colk = matC.getColumn (tCol);
      for (size_t i = 0, ie = colk.size (); i < ie; i++) {
        pppms[colk.keyAt (i)].setValue (tCol, 0);
      }
      colk.clear ();
      matB.getColumn (tCol).clear ();
    }

  private:
    static size_t sumAbsValues (const SparseArray<T> &col)
    {
      size_t tot = 0;
      for (size_t i = 0; i < col.size (); i++) {
        tot += std::abs (col.valueAt (i));
      }
      return tot;
    }

    static void test1b1 (MatrixCol<T> &matC, MatrixCol<T> &matB,
                         std::vector<PpPm> &pppms,
                         const Check11bResult &chkResult)
    {
      if (DEBUG) {
        std::cout << "Rule 1b.1 : " << chkResult.row << std::endl;
      }
      int tCol = chkResult.col;
      // [1.1.b.1] let k be the unique index of column belonging to P+ (resp. to P-)
      while (chkResult.p->size () > 0) {
        int j = chkResult.p->keyAt (0);
        // substitute to the column of index j the linear combination of
        // the columns indexed by k and j with the coefficients
        // |chj| and |chk| respectively.
        T chk = std::abs (matC.get (chkResult.row, tCol));
        T chj = std::abs (matC.get (chkResult.row, j));
        T gcdt = std::gcd (chk, chj);
        chk /= gcdt;
        chj /= gcdt;

        SparseBoolArray changed = sumProdInto (chk, matC.getColumn (j), chj,
                                               matC.getColumn (tCol));
        for (size_t ind = 0, inde = changed.size (); ind < inde; ind++) {
          pppms[changed.keyAt (ind)].setValue (
              j, matC.getColumn (j).get (changed.keyAt (ind)));
        }
        SparseArray<T> &coljb = matB.getColumn (j);
        sumProdInto (chk, coljb, chj, matB.getColumn (tCol));

        T gcdm = gcd(matC.getColumn (j));
        if (gcdm != 1) {
          T gcdb = gcd(matB.getColumn (j));
          if (gcdb != 1) {
            T ggcd = std::gcd(gcdm,gcdb);
            if (ggcd != 1) {
              matC.getColumn(j).scalarDiv(ggcd);
              matB.getColumn (j).scalarDiv(ggcd);
              //std::cout << "reduce by " << ggcd << std::endl;
            }
          }
        }
      }
      // delete from the extended matrix the column of index k
      clearColumn (chkResult.col, matC, matB, pppms);
    }

    static T gcd (const std::vector<T> & set)
    {
      if (set.size () == 0) return 0;
      T gcd = set[0];
      for (size_t i = 1; i < set.size (); i++) {
        gcd = InvariantCalculator::gcd (gcd, set[i]);
        if (gcd == 1) return 1;
      }
      return gcd;
    }

    static T gcd (const SparseArray<T> & set)
    {
      if (set.size () == 0) return 0;
      T gcd = set.valueAt (0);
      for (size_t i = 1; i < set.size (); i++) {
        gcd = std::gcd (gcd, set.valueAt (i));
        if (gcd == 1) return 1;
      }
      return gcd;
    }

    static void normalize (std::vector<T> & invariants)
    {
      int gcd = InvariantCalculator::gcd (invariants);
      if (gcd > 1) {
        for (size_t j = 0; j < invariants.size (); ++j) {
          invariants[j] /= gcd;
        }
      }
    }

  public:
    static void normalizeWithSign (SparseArray<T> &col)
    {
      bool allneg = true;
      for (size_t i = 0; i < col.size (); i++) {
        if (col.valueAt (i) > 0) {
          allneg = false;
          break;
        }
      }
      if (allneg) {
        for (size_t i = 0; i < col.size (); i++) {
          col.setValueAt (i, -col.valueAt (i));
        }
      }

      T gcd = InvariantCalculator::gcd (col);
      if (gcd > 1) {
        for (size_t j = 0; j < col.size (); ++j) {
          T norm = col.valueAt (j) / gcd;
          col.setValueAt (j, norm);
        }
        if (DEBUG)
          std::cout << "Applied gcd to invariant with gcd =" << gcd << std::endl;
      }
    }

    static void normalize (SparseArray<T> & invariants)
    {
      T gcd = InvariantCalculator::gcd (invariants);
      if (gcd > 1) {
        for (size_t j = 0; j < invariants.size (); ++j) {
          T norm = invariants.valueAt (j) / gcd;
          invariants.setValueAt (j, norm);
        }
      }
    }


  };

#endif /* INVARIANTCALCULATOR_H_ */
