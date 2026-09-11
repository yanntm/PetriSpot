#pragma once
#include <unordered_map>
#include "reduction/TransitionAlgebra.h"
namespace petri::reduction {
/** Transcription of Java ruleRedundantCompositionsBounds, outside reduce().
 * Group equal positive effects and discard stronger presets. This method is
 * NOT scheduled: its sole UpperBoundsSolver call is commented out in Java.
 * Its negative-effect omissions require a separate soundness assessment before
 * enabling it for arbitrary bound forms; preserving reference behavior means
 * retaining the inactive status rather than silently adding it to reachability. */
struct BoundsDominance {
  static constexpr const char* name = "bounds-dominance";
  template<class T> static void apply(Workspace<T>& w) {
    if (static_cast<size_t>(std::count(w.liveT.begin(), w.liveT.end(), true)) > w.config.redundantTransitionLimit) return;
    struct Bucket { SparseArray<T> effect; std::vector<size_t> ids; };
    std::unordered_map<size_t,std::vector<Bucket>> effects;
    std::vector<bool> removed(w.transitions.size());
    for (size_t t = 0; t < w.transitions.size(); ++t) if (w.liveT[t]) {
      if (w.stop()) return;
      auto effect = transitionEffect(w,t); SparseArray<T> positive;
      for (size_t i = 0; i < effect.size(); ++i) if (effect.valueAt(i) > 0) positive.put(effect.keyAt(i),effect.valueAt(i));
      auto& buckets = effects[positive.hash()];
      auto it = std::find_if(buckets.begin(), buckets.end(), [&](const auto& b) { return b.effect == positive; });
      if (it == buckets.end()) { buckets.push_back({std::move(positive),{t}}); continue; }
      for (size_t other : it->ids) {
        if (covers(w.pre.getColumn(t),w.pre.getColumn(other))) removed[t] = true;
        else if (covers(w.pre.getColumn(other),w.pre.getColumn(t))) removed[other] = true;
      }
      it->ids.push_back(t); std::erase_if(it->ids,[&](size_t t) { return removed[t]; });
    }
    for (size_t t = 0; t < removed.size(); ++t) if (removed[t]) w.retireTransition(t);
  }
};
}
