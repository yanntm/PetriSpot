// Tiny matrix-level regression: one named matrix per invocation, no net engine.
#include <cassert>
#include <iostream>
#include <random>
#include <string>
#include "invariants/InvariantMiddle.h"

using Matrix = MatrixCol<long>;
using Middle = petri::InvariantMiddle<long>;

static void verify(const Matrix& input, const Matrix& vectors, int direction) {
    assert(vectors.getRowCount() == input.getRowCount());
    for (const auto& b : vectors.getColumns()) {
        assert(b.size() != 0);
        for (size_t i = 0; i < b.size(); ++i) assert(b.valueAt(i) > 0);
        bool strict = false;
        for (const auto& constraint : input.getColumns()) {
            long effect = 0;
            for (size_t i = 0; i < constraint.size(); ++i)
                effect += constraint.valueAt(i) * b.get(constraint.keyAt(i));
            assert(direction * effect >= 0);
            strict |= effect != 0;
        }
        assert(strict);
    }
}

static Matrix example(const std::string& name) {
    if (name == "empty") return Matrix(4, 0);
    if (name == "sieve" || name == "growing") {
        Matrix m(4, 3);
        long sign = name == "sieve" ? -1 : 1;
        m.getColumn(0).put(2, sign);
        m.getColumn(1).put(3, sign);
        m.getColumn(2).put(3, sign);
        return m;
    }
    if (name == "mixed") {
        Matrix m(1, 2);
        m.getColumn(0).put(0, 1);
        m.getColumn(1).put(0, -1);
        return m;
    }
    if (name == "scaled") {
        Matrix m(3, 4);
        m.getColumn(0).put(0, -2);
        m.getColumn(0).put(1, -4);
        m.getColumn(1).put(0, -3);
        m.getColumn(1).put(1, -6);
        m.getColumn(2).put(2, 5);
        return m;
    }
    std::mt19937 random(static_cast<unsigned>(std::stoul(name)));
    Matrix m(6, 4);
    for (size_t j = 0; j < 4; ++j)
        for (size_t i = 0; i < 6; ++i) {
            long value = static_cast<long>(random() % 5) - 2;
            if (value) m.getColumn(j).put(i, value);
        }
    return m;
}

int main(int argc, char** argv) {
    assert(argc == 2);
    std::string name = argv[1];
    Matrix input = example(name);
    for (bool transpose : {false, true}) {
        Matrix m = transpose ? input.transpose() : input;
        for (bool semi : {false, true}) {
            for (int variant = 0; variant < 4; ++variant) {
                petri::EliminationHeuristic heur(variant != 1,
                    petri::EliminationHeuristic::PivotStrategy::FindBest,
                    -1, variant != 2, variant == 3, variant == 3, variant == 3);
                auto legacy = Middle::computePInvariants(m, semi, heur);
                auto collected = Middle::computePInvariantsWithInequalities(m, semi, heur);
                assert(legacy.first == collected.basis);
                assert(legacy.second.size() == collected.permutations.size());
                for (size_t i = 0; i < legacy.second.size(); ++i) {
                    assert(legacy.second[i].index == collected.permutations[i].index);
                    assert(legacy.second[i].elements == collected.permutations[i].elements);
                }
                verify(m, collected.inequalities.decreasing, -1);
                verify(m, collected.inequalities.increasing, 1);
                if (!transpose && name == "sieve") {
                    assert(collected.inequalities.decreasing.getColumnCount() == 2);
                    assert(collected.inequalities.increasing.getColumnCount() == 0);
                }
                if (!transpose && name == "growing") {
                    assert(collected.inequalities.increasing.getColumnCount() == 2);
                    assert(collected.inequalities.decreasing.getColumnCount() == 0);
                }
                if (!transpose && name == "mixed") {
                    assert(collected.inequalities.increasing.getColumnCount() == 0);
                    assert(collected.inequalities.decreasing.getColumnCount() == 0);
                }
            }
        }
    }
    std::cout << "PASS " << name << '\n';
}
