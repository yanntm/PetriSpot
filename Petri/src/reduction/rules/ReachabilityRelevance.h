#pragma once
#include "reduction/Workspace.h"
namespace petri::reduction {
/** Backward closure of observed coordinates: every transition changing a kept
 * place is kept, and so are all its input places. This stronger-than-necessary
 * dependency graph includes consuming effects. Projected paths lift because
 * all guards of a kept transition are retained; omitted steps stutter on the
 * closure. REACHABILITY only, with an explicit observation set (possibly empty).
 * O(P+T+arcs); no place cross-product graph is materialized. */
struct ReachabilityRelevance {
  static constexpr const char* name = "reachability-relevance";
  template<class T> static void apply(Workspace<T>& w) {
    if (w.config.goal != Goal::REACHABILITY || !w.config.relevance) return;
    std::vector<bool> keepP = w.observed, keepT(w.transitions.size(), false);
    std::vector<size_t> todo;
    for (size_t p = 0; p < keepP.size(); ++p) if (keepP[p]) todo.push_back(p);
    for (size_t k = 0; k < todo.size(); ++k) {
      if (k % 256 == 0 && w.stop()) return;
      size_t p = todo[k];
      for (const auto* matrix : {&w.consumers, &w.producers}) {
        const auto& row = matrix->getColumn(p);
        for (size_t j = 0; j < row.size(); ++j) {
          size_t t = row.keyAt(j);
          if (keepT[t] || w.pre.getColumn(t).get(p) == w.post.getColumn(t).get(p)) continue;
          keepT[t] = true;
          const auto& pre = w.pre.getColumn(t);
          for (size_t i = 0; i < pre.size(); ++i) {
            size_t q = pre.keyAt(i);
            if (!keepP[q]) { keepP[q] = true; todo.push_back(q); }
          }
        }
      }
    }
    for (size_t t = 0; t < keepT.size(); ++t) if (w.liveT[t] && !keepT[t]) w.retireTransition(t);
    for (size_t p = 0; p < keepP.size(); ++p) if (w.liveP[p] && !keepP[p]) w.retirePlace(p);
  }
};
}
