/*
 * MatrixCol.h
 *
 *  Created on: Oct 13, 2020
 *      Author: ythierry
 */

#ifndef MATRIXCOL_H_
#define MATRIXCOL_H_

/*-
 * APT - Analysis of Petri Nets and labeled Transition systems
 * Copyright (C) 2012-2013  Members of the project group APT
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General  License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General  License for more details.
 *
 * You should have received a copy of the GNU General  License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

// package fr.lip6.move.gal.util;
#include "SparseArray.h"
#include "InvariantHelpers.h"
#include <unordered_set>

/**
 * A Matrix specified for the invariant module, stored by COLUMN so that deletColumn has good complexity.
 * @author Manuel Gieseking, Dennis-Michael Borde, Yann Thierry-Mieg
 */
template<typename T>
  class MatrixCol
  {

    size_t iRows;
    size_t iCols;
    std::vector<SparseArray<T>> lCols;
  public:
    /**
     * Returns the identity matrix with given col and row count.
     * @param rows - the row count of the identity matrix.
     * @param cols - the col count of the identity matrix.
     * @return the identity matrix with the given col and row count.
     */
    static MatrixCol identity (size_t rows, size_t cols)
    {
      MatrixCol result (rows, cols);
      for (size_t i = 0; i < rows && i < cols; ++i) {
        result.set (i, i, 1);
      }
      return result;
    }

    /**
     * Constructor for a new Matrix with the given col and row count.
     * @param rows - the row count of the resulting matrix.
     * @param cols - the col count of the resulting matrix.
     */
    MatrixCol (size_t rows = 0, size_t cols = 0)
    {
      this->iRows = rows;
      this->iCols = cols;
      this->lCols.reserve (this->iCols);

      for (size_t col = 0; col < iCols; col++) {
        lCols.emplace_back ();
      }
    }

    /**
     * build a copy of a MatrixCol
     * @param ori
     */
    MatrixCol (const MatrixCol &ori)
    {
      iRows = ori.iRows;
      iCols = ori.iCols;
      lCols.reserve (iCols);
      for (const SparseArray<T> &a : ori.lCols) {
        lCols.emplace_back (a);
      }
    }

    MatrixCol& operator= (const MatrixCol &source)
    {
      if (this != &source) {
        lCols.clear ();
        iRows = source.iRows;
        iCols = source.iCols;
        lCols.reserve (iCols);
        for (const SparseArray<T> &a : source.lCols) {
          lCols.emplace_back (a);
        }
      }
      return *this;
    }

    MatrixCol (MatrixCol &&other) noexcept
        : iRows (other.iRows), iCols (other.iCols), lCols (
            std::move (other.lCols))
    {
      other.iRows = 0;
      other.iCols = 0;
    }

    MatrixCol& operator= (MatrixCol &&other) noexcept
    {
      if (this != &other) {
        iRows = other.iRows;
        iCols = other.iCols;
        lCols = std::move (other.lCols);
        other.iRows = 0;
        other.iCols = 0;
      }
      return *this;
    }

    /**
     * Constructor for a new Matrix with the values from the given array.
     * @param src - the template to create the matrix from.
     */
    MatrixCol (const int *const src, size_t rows, size_t cols)
    {
      this->iRows = rows;
      this->iCols = cols;
      this->lCols.reserve (this->iCols);

      for (size_t col = 0; col < this->iCols; ++col) {
        lCols.push_back (SparseArray<T> ());
        SparseArray<T> &toadd = lCols.back ();
        for (size_t row = 0; row < this->iRows; ++row) {
          int val = src[row + col * iRows];
          if (val != 0) {
            toadd.append (row, val);
          }
        }
      }
    }

    // Set columns using an rvalue reference.
    void setColumns (std::vector<SparseArray<T>> &&newCols)
    {
      lCols = std::move (newCols);
      iCols = lCols.size ();
    }

    void reserveColumns (size_t nbCol)
    {
      this->lCols.reserve (nbCol);
    }

    /** Exchange back to explicit form if required.
     * @return an explicit verion of this matrix, mat[i] giving a column at index i.
     */
    int* toExplicit () const
    {
      int *mat = new int[iRows * iCols];
      ::memset (mat, 0, iRows * iCols * sizeof(int));
      for (size_t col = 0; col < this->iCols; ++col) {
        const SparseArray<T> &arr = lCols[col];
        for (size_t i = 0; i < arr.size (); i++) {
          size_t row = arr.keyAt (i);
          int val = arr.valueAt (i);
          mat[row + col * iRows] = val;
        }
      }
      return mat;
    }

    /**
     * Returns the count of rows of this matrix.
     * @return the count of rows of this matrix.
     */
    size_t getRowCount () const
    {
      return this->iRows;
    }

    /**
     * Add a new empty row.
     */
    void addRow ()
    {
      iRows++;
    }

    /**
     * Returns the count of cols of this matrix.
     * @return the count of cols of this matrix.
     */
    size_t getColumnCount () const
    {
      return this->iCols;
    }

    size_t getEntryCount ()
    {
      size_t nbArcs = 0;
      for (const auto &col : getColumns ()) {
        nbArcs += col.size ();
      }
      return nbArcs;
    }

    /**
     * Returns the column with the given index of this matrix.
     * @param i - the index of the wished column of this matrix.
     * @return a copy of the column with the given index of this matrix.
     */
    SparseArray<T>& getColumn (size_t i)
    {
      return lCols[i];
    }
    const SparseArray<T>& getColumn (size_t i) const
    {
      return lCols[i];
    }

    SparseArray<T>& setColumn (size_t i, const SparseArray<T> &v)
    {
      return lCols[i] = std::move (v);
    }

    T get (size_t row, size_t col) const
    {
      return lCols[col].get (row);
    }

    void set (size_t row, size_t col, T val)
    {
      assert(row < iRows && col < iCols);
      if (val != 0) {
        SparseArray<T> &column = lCols[col];
        if (column.size () == 0 || column.keyAt (column.size () - 1) < row) {
          column.append (row, val);
        } else {
          column.put (row, val);
        }
      } else {
        lCols[col].del (row);
      }
    }

    /**
     * Deletes the column with the given index from this matrix.
     * @param j - the index of the column which should be deleted.
     */
    void deleteColumn (size_t j)
    {
      lCols.erase (lCols.begin () + j);
      this->iCols -= 1;
    }

    /**
     * Appends a given column to this matrix. That means adding the given column from the right side to this matrix.
     * @param column - the column to append.
     */
    void appendColumn (const SparseArray<T> &column)
    {
      assert(column.size () == 0 || iRows > column.keyAt (column.size () - 1));
      lCols.push_back (column);
      this->iCols++;
    }

    void appendColumn (SparseArray<T> &&column)
    {
      assert(column.size () == 0 || iRows > column.keyAt (column.size () - 1));
      lCols.emplace_back (column);
      this->iCols++;
    }

    /**
     * Checks if this matrix only contains of components equal to zero.
     * @return true if this matrix has just components equal to zero.
     */
    bool isZero () const
    {
      for (const SparseArray<T> &row : this->lCols) {
        if (row.size () != 0) {
          return false;
        }
      }
      return true;
    }

    /**
     * Transpose the Matrix in a new copy.
     */
    MatrixCol transpose () const
    {
      MatrixCol tr (iCols, iRows);
      transposeTo (tr, false);
      return tr;
    }

    T maxVal () const
    {
      if (lCols.empty ()) {
        return 0;
      } else {
        T max = 0;
        for (const SparseArray<T> &col : lCols) {
          auto mv = col.maxVal ();
          if (std::abs (mv) > std::abs (max)) {
            max = mv;
          }
        }
        return max;
      }
    }

    void transposeTo (MatrixCol &tr) const
    {
      transposeTo (tr, true);
    }
    void transposeTo (MatrixCol &tr, bool clear) const
    {
      if (clear) tr.clear (getColumnCount (), getRowCount ());
      for (size_t tcol = 0; tcol < iCols; tcol++) {
        const SparseArray<T> &col = lCols[tcol];
        for (size_t k = 0; k < col.size (); k++) {
          size_t trow = col.keyAt (k);
          T val = col.valueAt (k);
          tr.set (tcol, trow, val);
        }
      }
    }

    /**
     * Normalizes and reduces this matrix in-place, keeping only unique, non-zero columns.
     * Columns are normalized with sign, empty columns are removed, and duplicates are eliminated.
     * Modifies this matrix directly, reducing iCols and lCols accordingly.
     */
    void normalizeAndReduce (bool withSign = false)
    {
      std::unordered_set<SparseArray<T>*> seen (iCols);
      size_t writePos = 0;

      for (size_t i = 0; i < iCols; ++i) {
        SparseArray<T> &col = lCols[i];
        if (col.size () == 0) continue;

        if (withSign) petri::normalizeWithSign (col);
        else petri::normalize (col);

        if (seen.find (&col) == seen.end ()) { // Not found
          if (writePos != i) {
            lCols[writePos] = std::move (col);
          }
          seen.insert (&lCols[writePos]);  // Insert the new location
          writePos++;
        }
      }

      lCols.resize (writePos);
      iCols = writePos;
    }

    void print (std::ostream &os) const
    {
      os << "Matrix{lCols=";
      bool first = true;
      for (const auto &col : lCols) {
        if (first) first = false;
        else os << ", ";
        col.print (os);
      }
      os << '}';
    }
    friend std::ostream& operator<< (std::ostream &os, const MatrixCol &m)
    {
      m.print (os);
      return os;
    }

    /**
     * Use with care !! Basically, just consider this is readonly.
     * @return our very own storage !
     */
    const std::vector<SparseArray<T>>& getColumns () const
    {
      return lCols;
    }

    void sortByColumnSize (bool descending = false)
    {
      std::sort (
          lCols.begin (), lCols.end (),
          [descending] (const SparseArray<T> &a, const SparseArray<T> &b) {
            size_t aSize = a.size();
            size_t bSize = b.size();

            // First criterion: size of the columns
                        if (aSize != bSize) {
                          return descending ? aSize > bSize : aSize < bSize;
                        }

                        // If sizes are equal, check if columns are empty
              if (aSize == 0) {
                return false;  // Maintain stability for empty columns
              }

              // Second criterion: index of the first key
              size_t aFirstKey = a.keyAt(0);
              size_t bFirstKey = b.keyAt(0);
              return aFirstKey < bFirstKey;// Ascending tiebreaker
            });
    }

    void clear (size_t rowCount, size_t colCount)
    {
      lCols.resize (colCount);
      iCols = colCount;
      iRows = rowCount;
      for (SparseArray<T> &s : lCols) {
        s.clear ();
      }
    }

    void deleteRow (size_t row)
    {
      for (SparseArray<T> &col : lCols) {
        col.deleteAndShift (row);
      }
      iRows -= 1;
    }

    void deleteRows (const std::vector<unsigned int> &todel)
    {
//		if (getColumnCount() >= 1000)
//			getColumns().parallelStream().unordered().forEach(col -> col.deleteAndShift(todel));
//		else
//			getColumns().stream().unordered().forEach(col -> col.deleteAndShift(todel));

      for (SparseArray<T> &col : lCols) {
        col.deleteAndShift (todel);
      }
      iRows -= todel.size ();
    }

    bool operator== (const MatrixCol &other) const
    {
      if (this == &other) return true;
      if (iCols != other.iCols) return false;
      if (iRows != other.iRows) return false;
      return lCols == other.lCols;
    }

    static MatrixCol sumProd (int alpha, const MatrixCol &ta, int beta,
                              const MatrixCol &tb)
    {
      // Matrices should be homogeneous dimensions for sum-product operation."
      assert(
          ta.getColumnCount () == tb.getColumnCount () && ta.getRowCount () == tb.getRowCount ());

      MatrixCol mat (ta.getRowCount (), ta.getColumnCount ());
      for (size_t col = 0, cole = ta.getColumnCount (); col < cole; col++) {
        mat.setColumn (
            col,
            SparseArray<T>::sumProd (alpha, ta.getColumn (col), beta,
                                     tb.getColumn (col)));
      }
      return mat;
    }

  };

#endif /* MATRIXCOL_H_ */
