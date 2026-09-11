#pragma once
#include "core/Arithmetic.hpp"
#include "reduction/Workspace.h"
namespace petri::reduction {
/** Restricted pre-agglomeration: p is empty and unobserved; its unique feeder
 * h only produces one token in p. h has a nonempty preset, consumes no p and
 * is invisible. Every input place of h has h as its only consumer, ensuring
 * persistence without conflict. Consumers of p take one and do not return it.
 * Replace each f by h;f, then retire h,p. h can be delayed until f needs its
 * token: no competing transition can take h's reserved inputs. The nonempty
 * preset witnesses divergence-freedom of h alone. Reachability/deadlock only. */
struct PreAgglo {
  static constexpr const char* name = "pre-single-producer";
  template<class T> static void apply(Workspace<T>& w) {
    if (!permitsAgglomeration(w.config)) return;
    for (size_t p = 0; p < w.places.size(); ++p) {
      if (p % 256 == 0 && w.stop()) return;
      if (!w.liveP[p] || w.observed[p] || w.marks[p] != 0) continue;
      const auto& in = w.producers.getColumn(p); const auto& out = w.consumers.getColumn(p);
      if (in.size() != 1 || in.valueAt(0) != 1 || out.size() == 0) continue;
      size_t h = in.keyAt(0);
      const auto& preH = w.pre.getColumn(h);
      if (w.post.getColumn(h).size() != 1 || preH.size() == 0 || preH.get(p) || w.visible(h)) continue;
      bool ok = true;
      for (size_t i = 0; i < preH.size(); ++i)
        if (w.consumers.getColumn(preH.keyAt(i)).size() != 1) { ok = false; break; }
      size_t arcs = 0;
      for (size_t i = 0; ok && i < out.size(); ++i) {
        size_t f = out.keyAt(i);
        if (out.valueAt(i) != 1 || w.post.getColumn(f).get(p)) { ok = false; break; }
        size_t added = preH.size() + w.pre.getColumn(f).size();
        if (added > w.config.maxComposedArcs || arcs > w.config.maxComposedArcs - added) { ok = false; break; }
        arcs += added;
      }
      if (!ok) continue;
      std::vector<size_t> tails;
      std::vector<SparseArray<T>> columns;
      std::vector<std::string> names;
      for (size_t i = 0; i < out.size(); ++i) {
        size_t f = out.keyAt(i);
        SparseArray<T> col = w.pre.getColumn(f); col.put(p, 0);
        for (size_t j = 0; j < preH.size(); ++j) {
          size_t q = preH.keyAt(j);
          col.put(q, petri::addExact(col.get(q), preH.valueAt(j)));
        }
        tails.push_back(f); columns.push_back(std::move(col)); names.push_back(w.composedName(h, f));
      }
      for (size_t i = 0; i < tails.size(); ++i) {
        w.replacePre(tails[i], std::move(columns[i])); w.transitions[tails[i]] = std::move(names[i]);
      }
      w.retireTransition(h); w.retirePlace(p);
    }
  }
};
}
