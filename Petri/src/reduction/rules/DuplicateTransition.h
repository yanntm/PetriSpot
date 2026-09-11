#pragma once
#include <unordered_map>
#include "reduction/Workspace.h"

namespace petri::reduction {
/** Identical pre/post pairs induce the same steps and enabling predicates.
 * Hash only to select a bucket, then compare both columns. Stable last
 * representative retains its name and, under a counting record, the
 * multiplicity of every transition fused into it. */
struct DuplicateTransition {
  static constexpr const char* name = "duplicate-transition";
  template<class T> static void apply(Workspace<T>& w) {
    std::unordered_map<size_t, std::vector<size_t>> buckets;
    for (size_t t = w.transitions.size(); t-- > 0;) if (w.liveT[t]) {
      if (t % 256 == 0 && w.stop()) return;
      const auto& pre = w.pre.getColumn(t); const auto& post = w.post.getColumn(t);
      auto& bucket = buckets[pre.hash() * 31 + post.hash()];
      size_t survivor = w.transitions.size();
      for (size_t other : bucket)
        if (pre == w.pre.getColumn(other) && post == w.post.getColumn(other)) { survivor = other; break; }
      if (survivor < w.transitions.size()) w.fuseTransition(t, survivor); else bucket.push_back(t);
    }
  }
};
}
