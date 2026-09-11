#pragma once
#include "reduction/Composition.h"
namespace petri::reduction {
/** Java ruleFreeAgglo. Empty unobserved p with invisible unit feeders producing
 * only p and unit consumers, no transition on both sides. Reachability allows
 * delaying a feeder to its consumer without the persistence required by pre.
 * Simple phase: one feeder, at most one feeder input. Complex removes those
 * effort restrictions. Full H x F products replace both sides. */
struct FreeAgglo {
  static constexpr const char* name = "free-agglomeration";
  template<class T> static void apply(Workspace<T>& w, bool complex = false) {
    if (w.config.goal != Goal::REACHABILITY || !w.config.agglomeration) return;
    for (size_t p = 0; p < w.places.size(); ++p) {
      if (w.stop()) return;
      if (!w.liveP[p] || w.observed[p] || w.marks[p] != 0) continue;
      const auto& in = w.producers.getColumn(p); const auto& out = w.consumers.getColumn(p);
      if (!complex && in.size() != 1) continue;
      bool ok = true;
      std::vector<size_t> hs, fs;
      for (size_t i = 0; i < in.size(); ++i) {
        size_t h = in.keyAt(i);
        if (in.valueAt(i) != 1 || w.post.getColumn(h).size() != 1 || w.pre.getColumn(h).get(p)
            || (!complex && w.pre.getColumn(h).size() > 1) || w.visible(h)) { ok = false; break; }
        hs.push_back(h);
      }
      for (size_t i = 0; ok && i < out.size(); ++i) {
        size_t f = out.keyAt(i);
        if (out.valueAt(i) != 1 || w.post.getColumn(f).get(p)) { ok = false; break; }
        fs.push_back(f);
      }
      if (ok) agglomerate(w, p, hs, fs);
    }
  }
};
struct ComplexFreeAgglo {
  static constexpr const char* name = "free-agglomeration-complex";
  template<class T> static void apply(Workspace<T>& w) { FreeAgglo::apply(w, true); }
};
}
