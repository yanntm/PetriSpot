#pragma once
#include <algorithm>
#include "lp/DeadTransitions.h"
#include "reduction/Workspace.h"
namespace petri::reduction {
/** Transitions the state equation proves never enabled leave without a trace:
 * no state, no arc, no token was theirs, so every goal and the counting record
 * survive. LIVENESS keeps them, they are its answer. Budgeted by
 * `config.deadMs` over the whole reduction; a pass cut short resumes from
 * where it stopped on the next call (lp/DeadTransitions.h). */
struct DeadTransition {
  static constexpr const char* name = "dead-transition";
  template<class T> static void apply(Workspace<T>& w) {
    if (w.config.goal == Goal::LIVENESS || w.config.deadMs <= 0) return;
    long remaining = w.config.deadMs - w.dead.spentMs;
    if (remaining <= 0) { w.dead.limited = true; return; }
    auto start = std::chrono::steady_clock::now();
    auto deadline = std::min(w.deadline, start + std::chrono::milliseconds(remaining));
    lp::DeadReport report;
    auto dead = lp::deadTransitions(w.pre, w.post, w.marks, w.liveT, w.dead.cursor, deadline, w.config.deadPivots, report);
    w.retireDeadTransitions(std::move(dead));
    w.dead.found += report.dead; w.dead.tested += report.tested; w.dead.solves += report.solves;
    w.dead.pivots += report.pivots; w.dead.limited |= report.limited; ++w.dead.passes;
    w.dead.spentMs += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count();
  }
};
}
