#pragma once
#include "reduction/Workspace.h"

namespace petri::reduction {
/** Constant coordinates: equal pre/post rows, or initially zero with no
 * positive effect. Consumers requiring more than the constant are dead.
 * Erasing the remaining tests preserves every firing step. Observed coordinates
 * retain their marking, as do weighted fused components whose constant sum
 * represents several original markings. Others disappear, their tokens
 * recorded when a counting record is attached. LIVENESS skips this rule to retain
 * dead-transition obligations. */
struct ConstantPlace {
  static constexpr const char* name = "constant-place";
  template<class T> static void apply(Workspace<T>& w) {
    if (w.config.goal == Goal::LIVENESS) return; // dead transitions must remain as liveness obligations
    std::vector<size_t> dead, constants;
    for (size_t p = w.places.size(); p-- > 0;) {
      if ((p % 256 == 0 && w.stop()) || !w.liveP[p]) { if (w.limited) return; continue; }
      const auto& in = w.producers.getColumn(p);
      const auto& out = w.consumers.getColumn(p);
      bool constant = in == out;
      if (!constant && w.marks[p] == 0) {
        constant = true;
        for (size_t i = 0; i < in.size(); ++i)
          if (in.valueAt(i) > out.get(in.keyAt(i))) { constant = false; break; }
      }
      if (!constant) continue;
      for (size_t i = 0; i < out.size(); ++i)
        if (out.valueAt(i) > w.marks[p]) dead.push_back(out.keyAt(i));
      if (w.observed[p]) w.erasePlaceArcs(p);
      else constants.push_back(p);
    }
    w.retireConstantPlaces(std::move(constants));
    w.retireDeadTransitions(std::move(dead));
  }
};
}
