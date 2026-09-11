#pragma once
#include "reduction/Composition.h"
namespace petri::reduction {
/** Java rulePartialFreeAgglo: keep p and its unique unit consumer, but replace
 * individual invisible unit feeders producing only p by h;f. Other feeders
 * remain available. Candidate places are collected before the pass and retired
 * transitions are deferred until its end, as in the reference. */
struct PartialFreeAgglo {
  static constexpr const char* name = "partial-free-agglomeration";
  template<class T> static void apply(Workspace<T>& w) {
    if (w.config.goal != Goal::REACHABILITY || !w.config.agglomeration) return;
    std::vector<bool> candidates(w.places.size());
    for (size_t t = 0; t < w.transitions.size(); ++t) if (w.liveT[t]) {
      const auto& post = w.post.getColumn(t);
      if (post.size() == 1 && post.valueAt(0) == 1) {
        size_t p = post.keyAt(0);
        if (w.marks[p] == 0 && !w.observed[p] && !w.visible(t)) candidates[p] = true;
      }
    }
    std::vector<size_t> removed;
    size_t added = 0;
    for (size_t p = 0; p < candidates.size(); ++p) if (candidates[p]) {
      if (w.stop()) break;
      auto in = w.producers.getColumn(p), out = w.consumers.getColumn(p);
      if (out.size() > 1 || (out.size() == 1 && (out.valueAt(0) != 1 || in.get(out.keyAt(0))))) continue;
      for (size_t i = 0; i < in.size(); ++i) {
        size_t h = in.keyAt(i);
        if (w.visible(h) || w.post.getColumn(h).size() != 1 || in.valueAt(i) != 1) continue;
        for (size_t j = 0; j < out.size(); ++j) {
          size_t f = out.keyAt(j);
          auto pre = composeColumn(w.pre.getColumn(h), w.pre.getColumn(f), T{1}, p);
          auto post = composeColumn(w.post.getColumn(h), w.post.getColumn(f), T{1}, p);
          auto name = w.composedName(h, f);
          w.appendTransition(std::move(pre), std::move(post), std::move(name)); ++added;
        }
        removed.push_back(h);
      }
    }
    if (added) for (size_t t : removed) w.retireTransition(t);
  }
};
}
