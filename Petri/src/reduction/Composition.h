#pragma once
#include <unordered_map>
#include "core/Arithmetic.hpp"
#include "reduction/Workspace.h"

namespace petri::reduction {
/** Sparse weighted sum, omitting the intermediate coordinate. Agglomeration
 * guards establish interchangeability; this implements Java agglomerateAround's
 * summed presets, rather than a general sequential-composition precondition. */
template<class T>
SparseArray<T> composeColumn(const SparseArray<T>& h, const SparseArray<T>& f, T count, size_t omitted) {
  SparseArray<T> result = h;
  result.put(omitted, 0);
  for (size_t i = 0; i < f.size(); ++i) if (f.keyAt(i) != omitted) {
    size_t p = f.keyAt(i);
    result.put(p, petri::addExact(result.get(p), petri::multiplyExact(count, f.valueAt(i))));
  }
  return result;
}

/** Prepare the whole H x F product and hash-deduplicate it before any mutation.
 * Initial tokens, when present, follow the unique F continuation. No snapshots
 * or transition products are captured for tracing on the disabled path. */
template<class T>
bool agglomerate(Workspace<T>& w, size_t p, const std::vector<size_t>& hs, const std::vector<size_t>& fs) {
  struct Transition { SparseArray<T> pre, post; std::string name; };
  std::vector<Transition> product;
  std::unordered_map<size_t, std::vector<size_t>> buckets;
  size_t arcs = 0;
  for (size_t h : hs) for (size_t f : fs) {
    if (w.stop()) return false;
    T count = w.post.getColumn(h).get(p) / w.pre.getColumn(f).get(p);
    auto pre = composeColumn(w.pre.getColumn(h), w.pre.getColumn(f), count, p);
    auto post = composeColumn(w.post.getColumn(h), w.post.getColumn(f), count, p);
    auto& bucket = buckets[pre.hash() * 31 + post.hash()];
    bool duplicate = false;
    for (size_t i : bucket) if (product[i].pre == pre && product[i].post == post) { duplicate = true; break; }
    if (duplicate) continue;
    size_t added = pre.size() + post.size();
    if (added > w.config.maxComposedArcs || arcs > w.config.maxComposedArcs - added) return false;
    arcs += added;
    bucket.push_back(product.size());
    product.push_back({std::move(pre), std::move(post), w.composedName(h, f)});
  }
  // Sparse marking delta: multiplication replaces Java's repeated firing loop.
  SparseArray<T> marking;
  if (w.marks[p] != 0) {
    size_t f = fs.at(0);
    T count = w.marks[p] / w.pre.getColumn(f).get(p);
    const auto& post = w.post.getColumn(f);
    for (size_t i = 0; i < post.size(); ++i) {
      size_t q = post.keyAt(i);
      marking.put(q, petri::addExact(w.marks[q], petri::multiplyExact(count, post.valueAt(i))));
    }
  }
  for (size_t i = 0; i < marking.size(); ++i) w.marks[marking.keyAt(i)] = marking.valueAt(i);
  w.marks[p] = 0;
  for (size_t h : hs) w.retireTransition(h);
  for (size_t f : fs) w.retireTransition(f);
  for (auto& t : product) w.appendTransition(std::move(t.pre), std::move(t.post), std::move(t.name));
  w.retirePlace(p);
  return true;
}
}
