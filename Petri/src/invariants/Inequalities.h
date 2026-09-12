#ifndef PETRI_INVARIANTS_INEQUALITIES_H
#define PETRI_INVARIANTS_INEQUALITIES_H

#include <unordered_set>
#include "core/MatrixCol.h"
#include "core/InvariantHelpers.h"

namespace petri {

/** Original-variable certificates; no completeness claim and no constants. */
template<typename T>
struct Inequalities {
    MatrixCol<T> decreasing;
    MatrixCol<T> increasing;

    explicit Inequalities(size_t variables = 0)
        : decreasing(variables, 0), increasing(variables, 0) {}
};

/** Read-only observer of discarded phase-1 pairs C[j] = A B[j]. */
template<typename T>
class InequalityCollector {
    Inequalities<T>& output;
    std::unordered_set<SparseArray<T>> decreasingSeen;
    std::unordered_set<SparseArray<T>> increasingSeen;

    static int uniformSign(const SparseArray<T>& v) {
        int sign = 0;
        for (size_t i = 0; i < v.size(); ++i) {
            T value = v.valueAt(i);
            if (value == 0) continue;
            int next = value > 0 ? 1 : -1;
            if (sign != 0 && sign != next) return 0;
            sign = next;
        }
        return sign;
    }

public:
    explicit InequalityCollector(Inequalities<T>& result) : output(result) {}

    void observe(const SparseArray<T>& b, const SparseArray<T>& c) {
        int effectSign = uniformSign(c);
        if (effectSign == 0) return;
        int coefficientSign = uniformSign(b);
        if (coefficientSign == 0) return;
        SparseArray<T> candidate = b;
        // A collector-only orientation overflow must not abort elimination.
        try {
            if (coefficientSign < 0) {
                for (size_t i = 0; i < candidate.size(); ++i)
                    candidate.setValueAt(i, multiplyExact(T(-1), candidate.valueAt(i)));
            }
        } catch (const std::overflow_error&) {
            return;
        }
        normalize(candidate);
        bool decreasing = effectSign * coefficientSign < 0;
        auto& seen = decreasing ? decreasingSeen : increasingSeen;
        if (seen.insert(candidate).second)
            (decreasing ? output.decreasing : output.increasing).appendColumn(candidate);
    }
};

} // namespace petri
#endif
