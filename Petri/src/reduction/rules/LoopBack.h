#pragma once
#include "reduction/graph/Dependency.h"
namespace petri::reduction {
/** Java ruleReducePlaces loop-back branch (REACHABILITY only). An empty
 * unobserved place has a unique feeder h and a consumer f exactly inverse to h.
 * Recompute the property prefix without f. If every output of h is outside it,
 * prune the irrelevant region. One application returns to the coordinator. */
struct LoopBack {
  static constexpr const char* name = "loop-back";
  template<class T> static void apply(Workspace<T>& w) {
    if (w.config.goal != Goal::REACHABILITY || !w.config.relevance) return;
    for (size_t p = 0; p < w.places.size(); ++p) {
      if (w.stop()) return;
      if (!w.liveP[p] || w.observed[p] || w.marks[p] != 0 || w.producers.getColumn(p).size() != 1) continue;
      size_t h = w.producers.getColumn(p).keyAt(0);
      const auto& out = w.consumers.getColumn(p);
      for (size_t i = 0; i < out.size(); ++i) {
        size_t f = out.keyAt(i);
        if (w.pre.getColumn(h) != w.post.getColumn(f) || w.post.getColumn(h) != w.pre.getColumn(f)) continue;
        std::vector<bool> skipped(w.transitions.size()), kept(w.places.size()); skipped[f] = true;
        auto graph = dependencyGraph(w, skipped);
        observationSeeds(w, kept); graph.prefix(kept, [&] { return w.stop(); });
        if (w.stop()) return;
        bool irrelevant = true;
        const auto& post = w.post.getColumn(h);
        for (size_t j = 0; j < post.size(); ++j) if (kept[post.keyAt(j)]) { irrelevant = false; break; }
        if (irrelevant) { dropOutsidePrefix(w, kept); return; }
        break;
      }
    }
  }
};
}
