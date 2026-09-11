#pragma once
#include "reduction/Composition.h"
namespace petri::reduction {
/** Java rulePreAgglo: empty unobserved p, unit arcs, disjoint H/F; every H
 * produces only p, is invisible and divergent-free, and is strongly quasi-
 * persistent (each of its input places has no other consumer). Multiple H
 * are permitted. The cheap phase requires either H or F to be a singleton;
 * the complex phase admits their Cartesian product. */
struct PreAgglo {
  static constexpr const char* name = "pre-agglomeration";
  template<class T> static void apply(Workspace<T>& w, bool complex = false) {
    if (!permitsAgglomeration(w.config)) return;
    for (size_t p = 0; p < w.places.size(); ++p) {
      if (w.stop()) return;
      if (!w.liveP[p] || w.observed[p] || w.marks[p] != 0) continue;
      const auto& in = w.producers.getColumn(p); const auto& out = w.consumers.getColumn(p);
      if ((in.size() == 0) || (out.size() == 0) || (!complex && in.size() > 1 && out.size() > 1)) continue;
      bool ok = true;
      std::vector<size_t> hs, fs;
      for (size_t i = 0; i < in.size(); ++i) {
        size_t h = in.keyAt(i);
        const auto& pre = w.pre.getColumn(h);
        if (in.valueAt(i) != 1 || w.post.getColumn(h).size() != 1 || (pre.size() == 0)
            || pre.get(p) || w.visible(h)) { ok = false; break; }
        for (size_t j = 0; j < pre.size(); ++j)
          if (w.consumers.getColumn(pre.keyAt(j)).size() > 1) { ok = false; break; }
        if (!ok) break;
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
struct ComplexPreAgglo {
  static constexpr const char* name = "pre-agglomeration-complex";
  template<class T> static void apply(Workspace<T>& w) { PreAgglo::apply(w, true); }
};
}
