#ifndef INVARIANTSTRIVIAL_H_
#define INVARIANTSTRIVIAL_H_

#include <vector>
#include <unordered_set>
#include <iostream>
#include <limits>
#include "SparseArray.h"
#include "MatrixCol.h"
#include "SparseBoolArray.h"
#include "Arithmetic.hpp"

namespace petri {

//---------------------------------------------------------------------
// Efficiently removes trivial columns (i.e. empty columns in matC) from both
// matC and matB. For each trivial column, its corresponding invariant from
// matB is moved into trivialInv.
//
// @param matC       The incidence matrix.
// @param matB       The transformation matrix.
// @param trivialInv Vector accumulating trivial invariants.
//---------------------------------------------------------------------
template<typename T>
void cullConstantColumns(MatrixCol<T>& matC, MatrixCol<T>& matB,
                         std::vector<SparseArray<T>>& trivialInv)
{
    std::vector<size_t> trivialIndexes;
    size_t colCount = matC.getColumnCount();

    // First pass: accumulate indexes of trivial columns.
    for (size_t col = 0; col < colCount; ++col) {
        if (matC.getColumn(col).size() == 0) {
            trivialIndexes.push_back(col);
        }
    }

    // Only rebuild if we found any trivial columns.
    if (!trivialIndexes.empty()) {
        std::cout << "Culling trivial invariants: removed "
                  << trivialIndexes.size() << " columns." << std::endl;

        // Prepare new vectors for non-trivial columns.
        std::vector<SparseArray<T>> newColsC;
        std::vector<SparseArray<T>> newColsB;
        newColsC.reserve(colCount - trivialIndexes.size());
        newColsB.reserve(colCount - trivialIndexes.size());

        size_t trivialPos = 0; // Pointer into trivialIndexes (sorted in increasing order)
        for (size_t col = 0; col < colCount; ++col) {
            if (trivialPos < trivialIndexes.size() && col == trivialIndexes[trivialPos]) {
                // Found a trivial column: move its corresponding invariant from matB.
                trivialInv.push_back(std::move(matB.getColumn(col)));
                ++trivialPos;
                // Skip this column from matC.
            } else {
                // Non-trivial column: move from both matrices.
                newColsC.push_back(std::move(matC.getColumn(col)));
                newColsB.push_back(std::move(matB.getColumn(col)));
            }
        }
        // Replace the old columns with the new (culled) ones.
        matC.setColumns(std::move(newColsC));
        matB.setColumns(std::move(newColsB));
    }
}

//---------------------------------------------------------------------
// Efficiently removes duplicate columns in matC.
// For each duplicate column k (i.e. a column equal to a previously seen column r),
// we compute the difference diff = matB(k) - matB(r) using SparseArray::sumProd,
// add it to trivialInv, and record k as a duplicate.
// After processing all columns, if any duplicates were found, we rebuild matC and matB
// (skipping the duplicate columns).
//
// @param matC       The incidence matrix.
// @param matB       The transformation matrix.
// @param trivialInv Vector accumulating trivial invariants.
//---------------------------------------------------------------------
template<typename T>
void cullDuplicateColumns(MatrixCol<T>& matC, MatrixCol<T>& matB,
                          std::vector<SparseArray<T>>& trivialInv)
{
    size_t colCount = matC.getColumnCount();
    std::vector<size_t> duplicateIndexes;  // collect duplicate column indices
    std::unordered_map<SparseArray<T>, size_t> repMap;

    // First pass: detect duplicates.
    for (size_t col = 0; col < colCount; ++col) {
        const SparseArray<T>& columnC = matC.getColumn(col);
        auto it = repMap.find(columnC);
        if (it == repMap.end()) {
            repMap.insert({columnC, col});
        } else {
            // Duplicate column found: record its index.
            duplicateIndexes.push_back(col);
            // Compute diff = matB(col) - matB(rep) using sumProd.
            SparseArray<T> diff = SparseArray<T>::sumProd(1, matB.getColumn(col),
                                                          -1, matB.getColumn(it->second));
            trivialInv.push_back(std::move(diff));
        }
    }

    // Only rebuild if we found duplicates.
    if (!duplicateIndexes.empty()) {
        std::cout << "Culling duplicate invariants: removed "
                  << duplicateIndexes.size() << " columns." << std::endl;

        // duplicateIndexes are produced in increasing order since col increases.
        std::vector<SparseArray<T>> newColsC;
        std::vector<SparseArray<T>> newColsB;
        newColsC.reserve(colCount - duplicateIndexes.size());
        newColsB.reserve(colCount - duplicateIndexes.size());

        size_t dupPos = 0; // pointer into duplicateIndexes
        for (size_t col = 0; col < colCount; ++col) {
            if (dupPos < duplicateIndexes.size() && duplicateIndexes[dupPos] == col) {
                ++dupPos; // skip this duplicate column
            } else {
                newColsC.push_back(std::move(matC.getColumn(col)));
                newColsB.push_back(std::move(matB.getColumn(col)));
            }
        }

        // Replace the old columns with the new (culled) ones.
        matC.setColumns(std::move(newColsC));
        matB.setColumns(std::move(newColsB));
    }
}

} // namespace petri

#endif // INVARIANTSTRIVIAL_H_
