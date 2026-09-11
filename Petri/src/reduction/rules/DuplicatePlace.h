#pragma once
#include <unordered_map>
#include "reduction/Workspace.h"
namespace petri::reduction {
/** Java ensureUnique for places. Seed representatives with protected places,
 * then scan backwards. Identical pre/post rows imply constant marking
 * differences: a place with at least the representative's marking is redundant.
 * A smaller unprotected place replaces an unprotected representative only.
 * Deferred deletion keeps recognition based on unchanged sparse rows. */
struct DuplicatePlace {
  static constexpr const char* name = "duplicate-place";
  template<class T> static void apply(Workspace<T>& w) {
    if (w.config.goal == Goal::STATESPACE) return;
    std::unordered_map<size_t,std::vector<size_t>> buckets;
    auto representative = [&](size_t p) -> size_t& {
      const auto& in = w.producers.getColumn(p); const auto& out = w.consumers.getColumn(p);
      auto& bucket = buckets[in.hash()*31+out.hash()];
      for (size_t& q : bucket) if (in == w.producers.getColumn(q) && out == w.consumers.getColumn(q)) return q;
      bucket.push_back(p); return bucket.back();
    };
    for (size_t p = 0; p < w.places.size(); ++p) if (w.liveP[p] && w.observed[p]) representative(p) = p;
    std::vector<bool> removed(w.places.size());
    for (size_t p = w.places.size(); p-- > 0;) if (w.liveP[p] && !w.observed[p]) {
      if (w.stop()) return;
      size_t& q = representative(p);
      if (q == p) continue;
      if (w.marks[p] >= w.marks[q]) removed[p] = true;
      else if (!w.observed[q]) { removed[q] = true; q = p; }
    }
    for (size_t p = 0; p < removed.size(); ++p) if (removed[p]) w.retirePlace(p);
  }
};
}
