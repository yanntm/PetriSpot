#pragma once
#include "reduction/Workspace.h"
namespace petri::reduction {
/** Java ruleReduceTrans sink branch. An invisible pure consumer only weakens
 * the marking. Omitting it preserves reachable observations by monotonicity;
 * it must remain for deadlock questions. */
struct SinkTransition {
  static constexpr const char* name = "sink-transition";
  template<class T> static void apply(Workspace<T>& w) {
    if (w.config.goal != Goal::REACHABILITY) return;
    for (size_t t = 0; t < w.transitions.size(); ++t) {
      if (w.stop()) return;
      if (w.liveT[t] && w.post.getColumn(t).size() == 0 && !w.visible(t)) w.retireTransition(t);
    }
  }
};
}
