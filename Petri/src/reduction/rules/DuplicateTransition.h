#pragma once
#include <unordered_map>
#include "reduction/Workspace.h"

namespace petri::reduction {
/** Identical pre/post pairs induce the same steps and enabling predicates.
 * Hash only to select a bucket, then compare both columns. Stable first
 * representative retains its name. Arc multiplicities are not reconstructed
 * by this increment, so STATESPACE keeps duplicate transitions. */
struct DuplicateTransition {
  static constexpr const char* name = "duplicate-transition";
  template<class T> static void apply(Workspace<T>& w) {
    if (w.config.goal == Goal::STATESPACE) return;
    std::unordered_map<size_t, std::vector<size_t>> buckets;
    for (size_t t = 0; t < w.transitions.size(); ++t) if (w.liveT[t]) {
      if (t % 256 == 0 && w.stop()) return;
      const auto& pre = w.pre.getColumn(t); const auto& post = w.post.getColumn(t);
      auto& bucket = buckets[pre.hash() * 31 + post.hash()];
      bool duplicate = false;
      for (size_t other : bucket)
        if (pre == w.pre.getColumn(other) && post == w.post.getColumn(other)) { duplicate = true; break; }
      if (duplicate) w.retireTransition(t); else bucket.push_back(t);
    }
  }
};
}
