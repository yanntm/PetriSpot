#pragma once
#include <unordered_map>
#include <unordered_set>
#include "reduction/TransitionAlgebra.h"
namespace petri::reduction {
/** Java ruleRedundantCompositions. First maintain weakest presets per identical
 * effect. Then find guaranteed two-step sequences from each transition's minimal
 * enabling marking; remove transitions with that combined effect and stronger
 * presets. Candidate second steps consume a place produced by the first. */
struct RedundantComposition {
  static constexpr const char* name = "redundant-composition";
  template<class T> static void apply(Workspace<T>& w) {
    if (w.config.goal == Goal::LIVENESS || w.config.goal == Goal::STATESPACE
        || static_cast<size_t>(std::count(w.liveT.begin(), w.liveT.end(), true)) > w.config.redundantTransitionLimit) return;
    struct Bucket { SparseArray<T> effect; std::vector<size_t> transitions; };
    std::unordered_map<size_t, std::vector<Bucket>> effects;
    std::vector<size_t> ids;
    std::vector<bool> removed(w.transitions.size());
    auto bucketFor = [&](const SparseArray<T>& effect) -> std::vector<size_t>& {
      auto& buckets = effects[effect.hash()];
      for (auto& b : buckets) if (b.effect == effect) return b.transitions;
      buckets.push_back({effect,{}}); return buckets.back().transitions;
    };
    for (size_t t = 0; t < w.transitions.size(); ++t) if (w.liveT[t]) {
      if (w.stop()) return;
      ids.push_back(t);
      auto& bucket = bucketFor(transitionEffect(w,t));
      for (size_t other : bucket) {
        if (covers(w.pre.getColumn(t), w.pre.getColumn(other))) removed[t] = true;
        else if (covers(w.pre.getColumn(other), w.pre.getColumn(t))) removed[other] = true;
      }
      bucket.push_back(t);
      std::erase_if(bucket,[&](size_t t) { return removed[t]; });
    }
    if (w.config.goal != Goal::LTL && w.config.goal != Goal::LI_LTL) {
      std::stable_sort(ids.begin(), ids.end(), [&](size_t a, size_t b) {
        return w.pre.getColumn(a).size()+w.post.getColumn(a).size() > w.pre.getColumn(b).size()+w.post.getColumn(b).size();
      });
      for (size_t t : ids) if (!removed[t]) {
        if (w.stop()) return;
        const auto& post = w.post.getColumn(t);
        auto effect = transitionEffect(w,t);
        std::unordered_set<size_t> candidates;
        for (size_t i = 0; i < post.size(); ++i) {
          const auto& row = w.consumers.getColumn(post.keyAt(i));
          for (size_t j = 0; j < row.size(); ++j) candidates.insert(row.keyAt(j));
        }
        for (size_t f : candidates) {
          if (w.stop()) return;
          if (removed[f] || ((w.config.goal == Goal::SI_LTL || w.config.goal == Goal::SI_CTL) && w.visible(t) && w.visible(f))) continue;
          if (!covers(post, w.pre.getColumn(f))) continue;
          auto chain = addColumns(effect, transitionEffect(w,f));
          auto found = effects.find(chain.hash());
          if (found == effects.end()) continue;
          for (auto& bucket : found->second) if (bucket.effect == chain) {
            for (size_t other : bucket.transitions)
              if (other != t && other != f && covers(w.pre.getColumn(other), w.pre.getColumn(t))) removed[other] = true;
            std::erase_if(bucket.transitions,[&](size_t t) { return removed[t]; });
          }
        }
      }
    }
    for (size_t t = 0; t < removed.size(); ++t) if (removed[t]) w.retireTransition(t);
  }
};
}
