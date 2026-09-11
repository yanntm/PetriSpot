#pragma once
#include "core/Arithmetic.hpp"
#include "reduction/Workspace.h"
namespace petri::reduction {
/** F-continuation with one consumer f, controlled only by empty unobserved p.
 * Every feeder h produces a positive multiple k of f's consumption and never
 * consumes p; f does not feed p and is invisible. Normalize every h to h;f^k.
 * f cannot be disabled except by another f and can always finish; moving its
 * invisible production earlier preserves observed reachability and deadlock.
 * Prepare all weighted columns before mutation so overflow leaves this match
 * untouched. Reuse feeder slots, update transposes locally, then retire f,p.
 * This restricted variant avoids cross-products; broader Java cases are pending. */
struct PostAgglo {
  static constexpr const char* name = "post-single-consumer";
  template<class T> static void apply(Workspace<T>& w) {
    if (!permitsAgglomeration(w.config)) return;
    for (size_t p = 0; p < w.places.size(); ++p) {
      if (p % 256 == 0 && w.stop()) return;
      if (!w.liveP[p] || w.observed[p] || w.marks[p] != 0) continue;
      const auto& out = w.consumers.getColumn(p); const auto& in = w.producers.getColumn(p);
      if (out.size() != 1 || in.size() == 0) continue;
      size_t f = out.keyAt(0); T weight = out.valueAt(0);
      if (w.pre.getColumn(f).size() != 1 || w.post.getColumn(f).get(p) != 0 || w.visible(f)) continue;
      bool ok = true;
      size_t arcs = 0;
      for (size_t i = 0; i < in.size(); ++i) {
        size_t h = in.keyAt(i);
        if (w.pre.getColumn(h).get(p) != 0 || in.valueAt(i) % weight != 0) { ok = false; break; }
        size_t added = w.post.getColumn(h).size() + w.post.getColumn(f).size();
        if (added > w.config.maxComposedArcs || arcs > w.config.maxComposedArcs - added) { ok = false; break; }
        arcs += added;
      }
      if (!ok) continue;
      std::vector<size_t> feeders;
      std::vector<SparseArray<T>> columns;
      std::vector<std::string> names;
      for (size_t i = 0; i < in.size(); ++i) {
        size_t h = in.keyAt(i); T k = in.valueAt(i) / weight;
        SparseArray<T> col = w.post.getColumn(h); col.put(p, 0);
        const auto& tail = w.post.getColumn(f);
        for (size_t j = 0; j < tail.size(); ++j) {
          size_t q = tail.keyAt(j);
          col.put(q, petri::addExact(col.get(q), petri::multiplyExact(k, tail.valueAt(j))));
        }
        feeders.push_back(h); columns.push_back(std::move(col)); names.push_back(w.composedName(h, f));
      }
      for (size_t i = 0; i < feeders.size(); ++i) {
        w.replacePost(feeders[i], std::move(columns[i])); w.transitions[feeders[i]] = std::move(names[i]);
      }
      w.retireTransition(f); w.retirePlace(p);
    }
  }
};
}
