#ifndef EliminationHeuristic_H
#define EliminationHeuristic_H

#include <limits>
#include <cstddef> // for size_t
// Note: ssize_t is assumed available. On non-POSIX systems, use an appropriate signed type.

namespace petri {

class EliminationHeuristic {
public:
    enum class PivotStrategy { FindBest, FindWorst, FindFirst };

    // Constructor:
    // - useSingleSignRow: default true.
    // - pivotStrategy: default FindBest.
    // - loopLimit: default -1 interpreted as infinity (std::numeric_limits<size_t>::max()).
    EliminationHeuristic(bool useSingleSignRow = true,
           PivotStrategy pivotStrategy = PivotStrategy::FindBest,
           ssize_t loopLimit = -1,
           bool useCulling = true,
           bool minimize = false)
        : useSingleSignRow_(useSingleSignRow),
          pivotStrategy_(pivotStrategy),
          loopLimit_((loopLimit == -1) ? std::numeric_limits<size_t>::max()
                                      : static_cast<size_t>(loopLimit)),
                                        useCulling_(useCulling),
                                        minimize_(minimize)
    {}

    bool useSingleSignRow() const { return useSingleSignRow_; }
    PivotStrategy getPivotStrategy() const { return pivotStrategy_; }
    size_t getLoopLimit() const { return loopLimit_; }
    bool useCulling() const { return useCulling_; }
    bool useMinimization() const { return minimize_; }
private:
    const bool useSingleSignRow_;
    const PivotStrategy pivotStrategy_;
    const size_t loopLimit_;
    const bool useCulling_;
    const bool minimize_;
};

} // namespace petri

#endif // CONFIG_H
