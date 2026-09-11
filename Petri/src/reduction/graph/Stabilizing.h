#pragma once
#include "core/Arithmetic.hpp"
#include "reduction/Workspace.h"
namespace petri::reduction {
struct Stabilizing { std::vector<bool> places, transitions; };
/** Java computeStabilizing fixed point; signed sparse effects are calculated
 * with checked arithmetic. Stabilizing does not mean disabled or constant. */
template<class T> Stabilizing computeStabilizing(Workspace<T>& w) {
  Stabilizing stable{std::vector<bool>(w.places.size()), std::vector<bool>(w.transitions.size())};
  MatrixCol<T> effects(w.places.size(), 0);
  bool positive = false;
  for (size_t t = 0; t < w.transitions.size(); ++t) {
    if (w.stop()) return stable;
    SparseArray<T> effect;
    T sum = 0;
    if (w.liveT[t]) {
      const auto& pre = w.pre.getColumn(t); const auto& post = w.post.getColumn(t);
      for (size_t i = 0; i < pre.size(); ++i) effect.put(pre.keyAt(i), -pre.valueAt(i));
      for (size_t i = 0; i < post.size(); ++i)
        effect.put(post.keyAt(i), petri::addExact(effect.get(post.keyAt(i)), post.valueAt(i)));
      for (size_t i = 0; i < effect.size(); ++i) sum = petri::addExact(sum, effect.valueAt(i));
      positive |= sum > 0; stable.transitions[t] = sum < 0;
    }
    effects.appendColumn(std::move(effect));
  }
  if (positive) std::fill(stable.transitions.begin(), stable.transitions.end(), false);
  auto rows = effects.transpose();
  bool changed;
  do {
    changed = false;
    for (size_t p = 0; p < w.places.size(); ++p) if (w.liveP[p] && !stable.places[p]) {
      if (w.stop()) return stable;
      const auto& row = rows.getColumn(p);
      bool fed = false;
      for (size_t i = 0; i < row.size(); ++i)
        if (!stable.transitions[row.keyAt(i)] && row.valueAt(i) > 0) { fed = true; break; }
      if (!fed) { stable.places[p] = true; changed = true; }
    }
    for (size_t t = 0; t < w.transitions.size(); ++t) if (w.liveT[t] && !stable.transitions[t]) {
      if (w.stop()) return stable;
      const auto& col = effects.getColumn(t);
      for (size_t i = 0; i < col.size(); ++i) if (col.valueAt(i) < 0 && stable.places[col.keyAt(i)]) {
        stable.transitions[t] = true; changed = true; break;
      }
    }
  } while (changed);
  return stable;
}
}
