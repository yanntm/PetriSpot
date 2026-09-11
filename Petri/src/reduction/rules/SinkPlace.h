#pragma once
#include "reduction/Workspace.h"
namespace petri::reduction {
/** An unobserved place with no consuming arcs constrains no enabling test.
 * Projection preserves observed reachability. Token/state counting is outside
 * this rule's contract. */
struct SinkPlace {
  static constexpr const char* name = "sink-place";
  template<class T> static void apply(Workspace<T>& w) {
    if (w.config.goal == Goal::STATESPACE) return;
    for (size_t p = w.places.size(); p-- > 0;) {
      if (p % 256 == 0 && w.stop()) return;
      if (w.liveP[p] && !w.observed[p] && w.consumers.getColumn(p).size() == 0) w.retirePlace(p);
    }
  }
};
}
