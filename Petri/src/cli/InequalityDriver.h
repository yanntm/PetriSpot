#ifndef PETRI_CLI_INEQUALITYDRIVER_H
#define PETRI_CLI_INEQUALITYDRIVER_H

#include <filesystem>
#include <stdexcept>
#include "cli/Options.h"
#include "invariants/InvariantMiddle.h"
#include "io/SparseMatrixIO.h"

namespace petri::cli {

/** Validate only the new opt-in scenario. No restrictions on legacy requests. */
inline void validateInequalityOptions(const Options& o) {
    if (!o.inequalities()) return;
    if (!o.invariants()) throw std::invalid_argument("Inequality collection requires a flow/semiflow flag.");
    bool exporting = !o.decreasingKERSFile.empty() || !o.increasingKERSFile.empty();
    if (exporting && (o.pflows || o.psemiflows) && (o.tflows || o.tsemiflows))
        throw std::invalid_argument("Inequality export requires a single P/T orientation per invocation.");
    std::vector<std::filesystem::path> paths;
    for (const auto& name : {o.basisKERSFile, o.decreasingKERSFile, o.increasingKERSFile}) {
        if (name.empty()) continue;
        auto path = std::filesystem::weakly_canonical(name);
        for (const auto& previous : paths)
            if (path == previous)
                throw std::invalid_argument("Basis and inequality output paths must be distinct.");
        paths.push_back(path);
    }
}

template<typename T>
void outputInequalities(const Options& o, const Inequalities<T>& inequalities,
                        const std::vector<std::string>* names,
                        const std::vector<T>* initial) {
    if (!o.decreasingKERSFile.empty() && !SparseMatrixIO<T>::write(inequalities.decreasing, o.decreasingKERSFile))
        throw std::runtime_error("Could not write decreasing inequalities.");
    if (!o.increasingKERSFile.empty() && !SparseMatrixIO<T>::write(inequalities.increasing, o.increasingKERSFile))
        throw std::runtime_error("Could not write increasing inequalities.");
    std::cout << "Collected " << inequalities.decreasing.getColumnCount() << " decreasing and "
              << inequalities.increasing.getColumnCount() << " increasing inequalities (incomplete).\n";
    if (o.quiet || !o.basisKERSFile.empty() || !o.decreasingKERSFile.empty() || !o.increasingKERSFile.empty()) return;
    for (bool decreasing : {true, false}) {
        const auto& matrix = decreasing ? inequalities.decreasing : inequalities.increasing;
        for (const auto& b : matrix.getColumns()) {
            if (names) {
                std::stringstream expression;
                const std::vector<T> noInitial;
                try {
                    T constant = InvariantMiddle<T>::printEquation(b, initial ? *initial : noInitial,
                                                                  *names, expression);
                    if (initial)
                        std::cout << "ineq : " << expression.str() << (decreasing ? " <= " : " >= ")
                                  << constant << (decreasing ? " (nonincreasing)\n" : " (nondecreasing)\n");
                    else
                        std::cout << "rythm : " << expression.str() << "; effect "
                                  << (decreasing ? "<= 0" : ">= 0") << " (nonzero)\n";
                } catch (const std::overflow_error&) {
                    std::cerr << "Inequality constant overflow; coefficients remain available through API/KERS.\n";
                }
                continue;
            }
            std::cout << "ineq : M^T b " << (decreasing ? "<=" : ">=") << " 0; b = {";
            for (size_t i = 0; i < b.size(); ++i)
                std::cout << (i ? ", " : "") << "x" << b.keyAt(i) << ":" << b.valueAt(i);
            std::cout << "}\n";
        }
    }
}

/** One dispatch outside the solver; disabled requests use the original API. */
template<typename T>
typename InvariantMiddle<T>::Invariants computeRequestedInvariants(
    const Options& o, MatrixCol<T>& matrix, bool semi, const EliminationHeuristic& heur,
    const std::vector<std::string>* names = nullptr, const std::vector<T>* initial = nullptr) {
    if (!o.inequalities())
        return InvariantMiddle<T>::computePInvariants(matrix, semi, o.timeout, heur);
    auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(o.timeout);
    auto result = InvariantMiddle<T>::computePInvariantsWithInequalities(matrix, semi, heur, deadline);
    outputInequalities(o, result.inequalities, names, initial);
    return {std::move(result.basis), std::move(result.permutations)};
}

} // namespace petri::cli
#endif
