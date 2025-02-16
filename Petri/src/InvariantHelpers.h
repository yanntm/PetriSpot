#ifndef INVARIANTHELPERS_H_
#define INVARIANTHELPERS_H_

#include <numeric>
#include <stdexcept> // for std::overflow_error
#include <string>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <limits>
#include <cstdlib>
#include <cassert>
#include "SparseArray.h"
#include "MatrixCol.h"
#include "SparseBoolArray.h"
#include "Arithmetic.hpp"

namespace petri {

//---------------------------------------------------------------------
// Returns a SparseBoolArray of changed keys while computing
// the linear combination: ta = alpha * ta + beta * tb.
//---------------------------------------------------------------------
template<typename T>
SparseBoolArray sumProdInto(int alpha, SparseArray<T>& ta, int beta, const SparseArray<T>& tb)
{
    SparseBoolArray changed;
    SparseArray<T> flow(std::max(ta.size(), tb.size()));

    size_t i = 0;
    size_t j = 0;
    while (i < ta.size() || j < tb.size()) {
        unsigned int ki = (i == ta.size()) ?
            std::numeric_limits<unsigned int>::max() : ta.keyAt(i);
        unsigned int kj = (j == tb.size()) ?
            std::numeric_limits<unsigned int>::max() : tb.keyAt(j);
        if (ki == kj) {
            T val = petri::addExact(petri::multiplyExact(alpha, ta.valueAt(i)),
                                    petri::multiplyExact(beta, tb.valueAt(j)));
            if (val != 0) {
                flow.append(ki, val);
            }
            if (val != ta.valueAt(i)) {
                changed.set(ki);
            }
            i++;
            j++;
        } else if (ki < kj) {
            T val = petri::multiplyExact(alpha, ta.valueAt(i));
            if (val != 0) {
                flow.append(ki, val);
            }
            if (val != ta.valueAt(i)) {
                changed.set(ki);
            }
            i++;
        } else if (kj < ki) {
            T val = petri::multiplyExact(beta, tb.valueAt(j));
            if (val != 0) {
                flow.append(kj, val);
            }
            if (val != 0) {
                changed.set(kj);
            }
            j++;
        }
    }
    ta = std::move(flow);
    return changed;
}

//---------------------------------------------------------------------
// Returns the sign of the value.
//---------------------------------------------------------------------
template<typename T>
int signum(T value)
{
    return (value > 0) - (value < 0);
}

//---------------------------------------------------------------------
// Computes the greatest common divisor of all elements in a vector.
//---------------------------------------------------------------------
template<typename T>
T gcd(const std::vector<T>& set)
{
    if (set.size() == 0) return 0;
    T g = set[0];
    for (size_t i = 1; i < set.size(); i++) {
        g = std::gcd(g, set[i]);
        if (g == 1) return 1;
    }
    return g;
}

//---------------------------------------------------------------------
// Computes the greatest common divisor of all elements in a SparseArray.
//---------------------------------------------------------------------
template<typename T>
T gcd(const SparseArray<T>& set)
{
    if (set.size() == 0) return 0;
    T g = set.valueAt(0);
    for (size_t i = 1; i < set.size(); i++) {
        g = std::gcd(g, set.valueAt(i));
        if (g == 1) return 1;
    }
    return g;
}

//---------------------------------------------------------------------
// Normalizes a vector of invariants by dividing each element by the gcd.
//---------------------------------------------------------------------
template<typename T>
void normalize(std::vector<T>& invariants)
{
    T g = gcd(invariants);
    if (g > 1) {
        for (size_t j = 0; j < invariants.size(); ++j) {
            invariants[j] /= g;
        }
    }
}

//---------------------------------------------------------------------
// Normalizes a SparseArray with sign adjustment.
// If all entries are negative, their signs are flipped before normalization.
//---------------------------------------------------------------------
template<typename T>
void normalizeWithSign(SparseArray<T>& col)
{
    bool allneg = true;
    for (size_t i = 0; i < col.size(); i++) {
        if (col.valueAt(i) > 0) {
            allneg = false;
            break;
        }
    }
    if (allneg) {
        for (size_t i = 0; i < col.size(); i++) {
            col.setValueAt(i, -col.valueAt(i));
        }
    }

    T g = gcd(col);
    if (g > 1) {
        for (size_t j = 0; j < col.size(); ++j) {
            T norm = col.valueAt(j) / g;
            col.setValueAt(j, norm);
        }
        // Debug output commented out; enable if needed.
        // std::cout << "Applied gcd to invariant with gcd =" << g << std::endl;
    }
}

//---------------------------------------------------------------------
// Normalizes a SparseArray by dividing each element by the gcd.
//---------------------------------------------------------------------
template<typename T>
void normalize(SparseArray<T>& invariants)
{
    T g = gcd(invariants);
    if (g > 1) {
        for (size_t j = 0; j < invariants.size(); ++j) {
            T norm = invariants.valueAt(j) / g;
            invariants.setValueAt(j, norm);
        }
    }
}

//---------------------------------------------------------------------
// Computes the sum of the absolute values of the entries in a SparseArray.
//---------------------------------------------------------------------
template<typename T>
size_t sumAbsValues(const SparseArray<T>& col)
{
    size_t tot = 0;
    for (size_t i = 0; i < col.size(); i++) {
        tot += std::abs(col.valueAt(i));
    }
    return tot;
}

//---------------------------------------------------------------------
// Efficiently removes negative entries from a set of invariants.
// For each invariant in colsBsparse, if any value is negative, the invariant is erased.
//---------------------------------------------------------------------
template<typename T>
void removeNegativeValues(std::unordered_set<SparseArray<T>>& colsBsparse)
{
    for (auto it = colsBsparse.begin(); it != colsBsparse.end(); ) {
        const SparseArray<T>& a = *it;
        bool hasNegativeValue = false;
        for (size_t i = 0, ie = a.size(); i < ie; i++) {
            if (a.valueAt(i) < 0) {
                hasNegativeValue = true;
                break;
            }
        }
        if (hasNegativeValue) {
            it = colsBsparse.erase(it);
        } else {
            ++it;
        }
    }
}

} // namespace petri

#endif // INVARIANTHELPERS_H_
