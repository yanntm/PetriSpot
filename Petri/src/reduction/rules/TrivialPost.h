#pragma once
#include "reduction/Workspace.h"
namespace petri::reduction {
/** Invisible h;p;f chain, p initially empty, unit arcs, h only produces p,
 * f only consumes p, and exactly one producer/consumer at p. Replace h's
 * output by f's output and retire f. No column deletion, append or product.
 * Reachability and deadlock are preserved by expanding h to h;f; intermediate
 * states are unobserved. Temporal eligibility is deliberately deferred. */
struct TrivialPost {
  static constexpr const char* name = "trivial-post";
  template<class T> static void apply(Workspace<T>& w) {
    if (!permitsAgglomeration(w.config)) return;
    for (size_t p = 0; p < w.places.size(); ++p) {
      if (p % 256 == 0 && w.stop()) return;
      if (!w.liveP[p] || w.observed[p] || w.marks[p] != 0) continue;
      const auto& in = w.producers.getColumn(p); const auto& out = w.consumers.getColumn(p);
      if (in.size() != 1 || out.size() != 1 || in.valueAt(0) != 1 || out.valueAt(0) != 1) continue;
      size_t h = in.keyAt(0), f = out.keyAt(0);
      if (h == f || w.post.getColumn(h).size() != 1 || w.pre.getColumn(f).size() != 1
          || w.visible(h) || w.visible(f)) continue;
      auto result = w.post.getColumn(f);
      auto name = w.composedName(h, f);
      w.replacePost(h, std::move(result));
      w.transitions[h] = std::move(name);
      w.retireTransition(f); w.retirePlace(p);
    }
  }
};
}
