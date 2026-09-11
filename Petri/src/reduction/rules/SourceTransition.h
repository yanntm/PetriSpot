#pragma once
#include "reduction/Workspace.h"
namespace petri::reduction {
/** Java ruleReduceTrans source deduction: a transition with no input guard is
 * always enabled, proving deadlock freedom regardless of its outputs. */
struct SourceTransition {
  static constexpr const char* name = "source-transition";
  template<class T> static void apply(Workspace<T>& w) {
    if (w.config.goal != Goal::DEADLOCK) return;
    for (size_t t = 0; t < w.transitions.size(); ++t) {
      if (w.stop()) return;
      if (w.liveT[t] && w.pre.getColumn(t).size() == 0) { w.deadlock = false; return; }
    }
  }
};
}
