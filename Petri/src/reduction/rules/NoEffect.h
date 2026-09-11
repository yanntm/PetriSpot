#pragma once
#include "reduction/Workspace.h"
namespace petri::reduction {
/** Identity firings add no reachable valuations. They can prevent deadlocks
 * or affect next-step/temporal behavior, hence REACHABILITY, and STATESPACE
 * once the arcs are no longer tracked: while they are, a consumer still counts
 * the guards of these transitions, so they stay. */
struct NoEffect {
  static constexpr const char* name = "no-effect";
  template<class T> static void apply(Workspace<T>& w) {
    bool arcsTracked = w.counting && w.counting->tmult;
    if (w.config.goal != Goal::REACHABILITY && !(w.config.goal == Goal::STATESPACE && !arcsTracked)) return;
    for (size_t t = 0; t < w.transitions.size(); ++t) {
      if (t % 256 == 0 && w.stop()) return;
      if (w.liveT[t] && w.pre.getColumn(t) == w.post.getColumn(t)) w.retireTransition(t);
    }
  }
};
}
