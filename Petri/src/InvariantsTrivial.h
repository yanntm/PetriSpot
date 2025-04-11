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
size_t cullConstantColumns(MatrixCol<T>& matC, MatrixCol<T>& matB,
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

    size_t tocull = trivialIndexes.size();
    // Only rebuild if we found any trivial columns.
    if (tocull > 0) {
        std::cout << "Culling trivial invariants: removing "
                  << tocull << " columns." << std::endl;

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
    return tocull;
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
// Helper function to generate a unique string key for a SparseArray based on GCD normalization
template<typename T>
std::string gcdString(const SparseArray<T>& v) {
    if (v.size() == 0) {
        return "empty";
    }
    // Compute GCD of absolute values
    T gcd = std::abs(v.valueAt(0));
    for (size_t i = 1; i < v.size(); ++i) {
        gcd = std::gcd(gcd, std::abs(v.valueAt(i)));
        if (gcd == 1) break;
    }
    if (gcd == 0) {
        return "zero";
    }
    // Normalize sign based on the first non-zero value
    T first_val = v.valueAt(0);
    int sign = (first_val > 0) ? 1 : -1;
    // Create string representation of normalized values
    std::ostringstream oss;
    for (size_t i = 0; i < v.size(); ++i) {
        if (i > 0) oss << ",";
        T normalized_val = (v.valueAt(i) / gcd) * sign;
        oss << v.keyAt(i) << ":" << normalized_val;
    }
    return oss.str();
}

// Extended function to cull duplicate columns, including rational multiples
template<typename T>
size_t cullDuplicateColumns(MatrixCol<T>& matC, MatrixCol<T>& matB,
                            std::vector<SparseArray<T>>& trivialInv)
{
    size_t colCount = matC.getColumnCount();
    std::vector<size_t> duplicateIndexes;  // Collect duplicate column indices
    std::unordered_map<std::string, size_t> repMap;  // Map from gcdString to representative index

    // First pass: detect duplicates up to rational multiples
    for (size_t col = 0; col < colCount; ++col) {
        const SparseArray<T>& columnC = matC.getColumn(col);
        std::string key = gcdString(columnC);
        auto it = repMap.find(key);
        if (it == repMap.end()) {
            repMap[key] = col;
        } else {
            // Duplicate found (rational multiple)
            duplicateIndexes.push_back(col);
            size_t rep = it->second;
            const SparseArray<T>& repColumnC = matC.getColumn(rep);
            // Use first values as coefficients to compute difference
            T a = columnC.valueAt(0);
            T b = repColumnC.valueAt(0);
            // Compute diff = b * matB(col) - a * matB(rep)
            SparseArray<T> diff = petri::sumProd(b, matB.getColumn(col),
                                                          -a, matB.getColumn(rep));
            trivialInv.push_back(std::move(diff));
        }
    }

    size_t tocull = duplicateIndexes.size();
    // Rebuild matrices only if duplicates were found
    if (tocull > 0) {
        std::cout << "Culling duplicate invariants (modulo rational): removing "
                  << tocull << " columns." << std::endl;

        // Rebuild matC and matB, skipping duplicates
        std::vector<SparseArray<T>> newColsC;
        std::vector<SparseArray<T>> newColsB;
        newColsC.reserve(colCount - tocull);
        newColsB.reserve(colCount - tocull);

        size_t dupPos = 0; // Pointer into duplicateIndexes
        for (size_t col = 0; col < colCount; ++col) {
            if (dupPos < tocull && duplicateIndexes[dupPos] == col) {
                ++dupPos; // Skip this duplicate column
            } else {
                newColsC.push_back(std::move(matC.getColumn(col)));
                newColsB.push_back(std::move(matB.getColumn(col)));
            }
        }

        // Update matrices with culled columns
        matC.setColumns(std::move(newColsC));
        matB.setColumns(std::move(newColsB));
    }
    return tocull;
}


} // namespace petri

#endif // INVARIANTSTRIVIAL_H_
