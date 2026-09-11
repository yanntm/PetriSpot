#pragma once
#include "reduction/Workspace.h"
namespace petri::reduction {
/** Java testModuloIsomorphism: if both columns of b are k times those of a,
 * k>1, firing b can be expanded into k firings of a. Locate candidates through
 * input places with a non-unit arc. Pair orientation is local and immutable. */
struct ScalarTransition {
  static constexpr const char* name = "scalar-transition";
  template<class T> static bool multiple(const SparseArray<T>& a, const SparseArray<T>& b, T k) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i)
      if (a.keyAt(i) != b.keyAt(i) || b.valueAt(i) % a.valueAt(i) != 0 || b.valueAt(i) / a.valueAt(i) != k) return false;
    return true;
  }
  template<class T> static void apply(Workspace<T>& w) {
    if (w.config.goal == Goal::LIVENESS || w.config.goal == Goal::STATESPACE) return;
    std::vector<bool> removed(w.transitions.size());
    for (size_t p = 0; p < w.places.size(); ++p) if (w.liveP[p]) {
      if (w.stop()) return;
      const auto& row = w.consumers.getColumn(p);
      bool weighted = false;
      for (size_t i = 0; i < row.size(); ++i) weighted |= row.valueAt(i) > 1;
      if (!weighted) continue;
      for (size_t i = 0; i < row.size(); ++i) for (size_t j = i + 1; j < row.size(); ++j) {
        if (w.stop()) return;
        size_t a = row.keyAt(i), b = row.keyAt(j);
        T va = row.valueAt(i), vb = row.valueAt(j);
        if (va == vb) continue;
        if (va > vb) { std::swap(a,b); std::swap(va,vb); }
        if (vb % va) continue;
        T k = vb / va;
        if (multiple(w.pre.getColumn(a), w.pre.getColumn(b), k)
            && multiple(w.post.getColumn(a), w.post.getColumn(b), k)) removed[b] = true;
      }
    }
    for (size_t t = 0; t < removed.size(); ++t) if (removed[t]) w.retireTransition(t);
  }
};
}
