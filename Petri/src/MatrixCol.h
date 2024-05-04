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

#include "SparseIntArray.h"

/**
 * A Matrix specified for the invariant module, stored by COLUMN so that deletColumn has good complexity.
 * @author Manuel Gieseking, Dennis-Michael Borde, Yann Thierry-Mieg
 */
class MatrixCol {

	size_t iRows;
	size_t iCols;
	std::vector<SparseIntArray> lCols;
public :
	/**
	 * Returns the identity matrix with given col and row count.
	 * @param rows - the row count of the identity matrix.
	 * @param cols - the col count of the identity matrix.
	 * @return the identity matrix with the given col and row count.
	 */
	static MatrixCol identity(size_t rows, size_t cols) {
		MatrixCol result(rows, cols);
		for (size_t i = 0; i < rows && i < cols; ++i) {
			result.set(i,i, 1);
		}
		return result;
	}

	/**
	 * Constructor for a new Matrix with the given col and row count.
	 * @param rows - the row count of the resulting matrix.
	 * @param cols - the col count of the resulting matrix.
	 */
	 MatrixCol(size_t rows=0, size_t cols=0) {
		this->iRows = rows;
		this->iCols = cols;
		this->lCols.reserve(this->iCols);

		for (size_t col = 0 ; col < iCols ; col++) {
			lCols.emplace_back();
		}
	}

	/**
	 * build a copy of a MatrixCol
	 * @param ori
	 */
	 MatrixCol (const MatrixCol & ori) {
		iRows = ori.iRows;
		iCols = ori.iCols;
		lCols.reserve(iCols);
		for (const SparseIntArray & a : ori.lCols) {
			lCols.emplace_back(a);
		}
	}

	/**
	 * Constructor for a new Matrix with the values from the given array.
	 * @param src - the template to create the matrix from.
	 */
	 MatrixCol(const int*const src, size_t rows, size_t cols) {
		this->iRows = rows;
		this->iCols = cols;
		this->lCols.reserve(this->iCols);

		for (size_t col = 0; col < this->iCols; ++col) {
			lCols.push_back(SparseIntArray());
			SparseIntArray & toadd = lCols.back();
			for (size_t row = 0; row < this->iRows; ++row) {
				int val = src[row + col * iRows];
				if (val != 0) {
					toadd.put(row,val);
				}
			}
		}
	}

	 void reserveColumns (size_t nbCol) {
		 this->lCols.reserve(nbCol);
	 }

	/** Exchange back to explicit form if required.
	 * @return an explicit verion of this matrix, mat[i] giving a column at index i.
	 */
	 int * toExplicit () const {
		int  * mat = new int[iRows*iCols];
		::memset(mat,0,iRows*iCols*sizeof(int));
		for (size_t col = 0; col < this->iCols; ++col) {
			const SparseIntArray & arr = lCols[col];
			for (size_t i = 0 ; i < arr.size() ; i++) {
				size_t row = arr.keyAt(i);
				int val = arr.valueAt(i);
				mat [row + col*iRows] =  val;
			}
		}
		return mat;
	}


	/**
	 * Returns the count of rows of this matrix.
	 * @return the count of rows of this matrix.
	 */
	 size_t getRowCount() const {
		return this->iRows;
	}

	/**
	 * Add a new empty row.
	 */
	 void addRow() {
		iRows ++;
	}

	/**
	 * Returns the count of cols of this matrix.
	 * @return the count of cols of this matrix.
	 */
	 size_t getColumnCount() const {
		return this->iCols;
	}

	/**
	 * Returns the column with the given index of this matrix.
	 * @param i - the index of the wished column of this matrix.
	 * @return a copy of the column with the given index of this matrix.
	 */
	 SparseIntArray & getColumn(size_t i) {
		return lCols[i];
	}
	const SparseIntArray & getColumn(size_t i) const {
		return lCols[i];
	}

	 SparseIntArray & setColumn(size_t i, const SparseIntArray & v) {
		return lCols[i] = v;
	}

	 int get (size_t row, size_t col) const {
		return lCols[col].get(row);
	}

	 void set (size_t row, size_t col, int val) {
		assert (! (row < 0 || col < 0 || row >= iRows || col >= iCols));
		if (val != 0) {
			SparseIntArray & column = lCols[col];
			if (column.size()== 0 || column.keyAt(column.size()-1) < row) {
				column.append(row, val);
			} else {
				column.put(row,val);
			}
		} else {
			lCols[col].del(row);
		}
	}

	/**
	 * Returns a row which has at least one component different from zero. It also returns the index of the column
	 * where a component not equal to zero was found. If such a row does not exists, then -1,-1.
	 * @return the index of the column with a none zero component and the addicted row or null if not existent.
	 */
	 std::pair<int,int> getNoneZeroRow() const {
		// optimize to prefer to return 1
		for (size_t tcol = 0; tcol < getColumnCount(); tcol++) {
			if (lCols.at(tcol).size()==0) {
				continue;
			} else {
				return {lCols[tcol].keyAt(0) , tcol};
			}
		}
		return {-1,-1};
	}



	/**
	 * Deletes the column with the given index from this matrix.
	 * @param j - the index of the column which should be deleted.
	 */
	 void deleteColumn(size_t j) {
		lCols.erase(lCols.begin()+j);
		this->iCols -= 1;
	}

	/**
	 * Appends a given column to this matrix. That means adding the given column from the right side to this matrix.
	 * @param column - the column to append.
	 */
	 void appendColumn(const SparseIntArray & column) {
		assert (column.size()==0 || iRows > column.keyAt(column.size()-1));
		lCols.push_back(column);
		this->iCols++;
	}

	/**
	 * Checks if this matrix only contains of components equal to zero.
	 * @return true if this matrix has just components equal to zero.
	 */
	 bool isZero() const {
		for (const SparseIntArray & row : this->lCols) {
			if (row.size() != 0) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Transpose the Matrix in a new copy.
	 */
	 MatrixCol transpose() const {
		MatrixCol tr(iCols, iRows);
		transposeTo(tr,false);
		return tr;
	}

	 void transposeTo(MatrixCol & tr) const {
		transposeTo(tr, true);
	}
	 void transposeTo(MatrixCol & tr, bool clear) const {
		if (clear)
			tr.clear(getColumnCount(),getRowCount());
		for (size_t tcol = 0; tcol < iCols; tcol++) {
			const SparseIntArray & col = lCols[tcol];
			for (size_t k =0 ; k < col.size() ; k++) {
				size_t trow = col.keyAt(k);
				int val = col.valueAt(k);
				tr.set(tcol, trow, val);
			}
		}
	}


	 void print (std::ostream & os) const {
		 os << "Matrix{lCols=" ;
		 bool first = true;
		 for (const auto & col : lCols) {
			 if (first)
				 first=false;
			 else
				 os << ", ";
			 col.print(os);
		 }
		 os << '}';
	 }

	/**
	 * Use with care !! Basically, just consider this is readonly.
	 * @return our very own storage !
	 */
	 const std::vector<SparseIntArray> & getColumns () {
		return lCols;
	}

	 void clear(size_t rowCount, size_t colCount) {
		 lCols.resize(colCount);
		iCols = colCount;
		iRows = rowCount;
		for (SparseIntArray & s : lCols) {
			s.clear();
		}
	}

	 void deleteRow(size_t row) {
		for (SparseIntArray & col : lCols) {
			col.deleteAndShift(row);
		}
		iRows -= 1;
	}

	 void deleteRows(const std::vector<unsigned int> & todel) {
//		if (getColumnCount() >= 1000)
//			getColumns().parallelStream().unordered().forEach(col -> col.deleteAndShift(todel));
//		else
//			getColumns().stream().unordered().forEach(col -> col.deleteAndShift(todel));

		for (SparseIntArray & col : lCols) {
			 col.deleteAndShift(todel);
		}
		iRows -= todel.size();
	}

	 bool operator==(const MatrixCol & other) {
		if (this == &other)
			return true;
		if (iCols != other.iCols)
			return false;
		if (iRows != other.iRows)
			return false;
		return lCols == other.lCols;
	}

	 static MatrixCol sumProd (int alpha, const MatrixCol & ta, int beta, const MatrixCol & tb) {
		 //			throw new IllegalArgumentException("Matrices should be homogeneous dimensions for sum-product operation.");
		assert (ta.getColumnCount() != tb.getColumnCount() || ta.getRowCount() != tb.getRowCount()) ;

		MatrixCol mat(ta.getRowCount(), ta.getColumnCount());
		for (size_t col=0,cole=ta.getColumnCount(); col < cole ; col++) {
			mat.setColumn(col, SparseIntArray::sumProd(alpha, ta.getColumn(col), beta, tb.getColumn(col)));
		}
		return mat;
	}

};



#endif /* MATRIXCOL_H_ */
