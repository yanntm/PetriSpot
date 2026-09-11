#pragma once
#include <map>
#include "reduction/TransitionAlgebra.h"
namespace petri::reduction {
/** Java ruleFusePlaceByFuture. Group unobserved places by outgoing degree and
 * greedily match their consumers modulo the two coordinates. Redirect incoming
 * arcs to the representative and sum markings; remove the other's consumers.
 * This is future equivalence, not free-SCC token redistribution. */
struct FutureEquivalent {
  static constexpr const char* name = "future-equivalent";
  /** Literal equalUptoPerm recognition: matching nonzero entries at either
   * exchanged coordinate is forbidden; the two exceptional entries must agree.
   * Keep the reference's merge-walk stopping condition and consumer match order. */
  template<class T> static bool equalUptoPerm(const SparseArray<T>& a, const SparseArray<T>& b, size_t pi, size_t pj) {
    if (a.size() != b.size()) return false;
    if (a.size() == 0) return true;
    T seenA = -1, seenB = -1;
    size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
      size_t ka = a.keyAt(i), kb = b.keyAt(j);
      T va = a.valueAt(i), vb = b.valueAt(j);
      if (ka == kb) {
        if (ka == pi || ka == pj || va != vb) return false;
        ++i; ++j;
      } else if (ka == pi && kb == pj) {
        if (va != vb) return false;
        seenA = va; seenB = vb; ++i; ++j;
      } else if (ka < kb) {
        if (ka != pi) return false;
        if (seenA == -1) seenA = va;
        if (seenB != -1 && seenB != va) return false;
        ++i;
      } else {
        if (kb != pj) return false;
        if (seenB == -1) seenB = vb;
        if (seenA != -1 && seenA != vb) return false;
        ++j;
      }
    }
    return seenA == seenB;
  }
  template<class T> static void apply(Workspace<T>& w) {
    if (w.config.goal == Goal::LIVENESS || w.config.goal == Goal::STATESPACE) return;
    std::map<size_t,std::vector<size_t>> buckets;
    for (size_t p = 0; p < w.places.size(); ++p) if (w.liveP[p] && !w.observed[p])
      buckets[w.consumers.getColumn(p).size()].push_back(p);
    std::map<size_t,size_t> fusion;
    for (const auto& [degree, ids] : buckets) {
      if (ids.size() >= w.config.futureBucketLimit) continue;
      for (size_t i = 0; i < ids.size(); ++i) {
        size_t pi = ids[i]; if (fusion.contains(pi)) continue;
        const auto& a = w.consumers.getColumn(pi);
        for (size_t j = i+1; j < ids.size(); ++j) {
          if (w.stop()) return;
          size_t pj = ids[j]; if (fusion.contains(pj)) continue;
          const auto& b = w.consumers.getColumn(pj);
          std::vector<bool> matched(b.size());
          bool equal = true;
          for (size_t ti = 0; ti < a.size(); ++ti) {
            size_t t = a.keyAt(ti); bool found = false;
            for (size_t tj = 0; tj < b.size(); ++tj) {
              size_t u = b.keyAt(tj);
              if (t == u) break; // reference deliberately rejects this match path
              if (matched[tj]) continue;
              if (equalUptoPerm(w.post.getColumn(t),w.post.getColumn(u),pi,pj)
                  && equalUptoPerm(w.pre.getColumn(t),w.pre.getColumn(u),pi,pj)) {
                found = true; matched[tj] = true; break;
              }
            }
            if (!found) { equal = false; break; }
          }
          if (equal) fusion[pj] = pi;
        }
      }
    }
    // Prepare changed markings and post columns before modifying the net.
    std::map<size_t,T> marking;
    std::map<size_t,SparseArray<T>> posts;
    std::vector<bool> removed(w.transitions.size());
    for (const auto& [pj, pi] : fusion) {
      if (w.stop()) return;
      auto [mark, inserted] = marking.try_emplace(pi,w.marks[pi]);
      mark->second = petri::addExact(mark->second,w.marks[pj]);
      const auto& producers = w.producers.getColumn(pj);
      for (size_t i = 0; i < producers.size(); ++i) {
        size_t t = producers.keyAt(i);
        auto [it, fresh] = posts.try_emplace(t,w.post.getColumn(t));
        T value = it->second.get(pj); it->second.put(pj,0);
        it->second.put(pi,petri::addExact(it->second.get(pi),value));
      }
      const auto& consumers = w.consumers.getColumn(pj);
      for (size_t i = 0; i < consumers.size(); ++i) removed[consumers.keyAt(i)] = true;
    }
    for (const auto& [p,m] : marking) w.marks[p] = m;
    for (auto& [t,col] : posts) w.replacePost(t,std::move(col));
    for (size_t t = 0; t < removed.size(); ++t) if (removed[t]) w.retireTransition(t);
    for (const auto& [pj,pi] : fusion) w.retirePlace(pj);
  }
};
}
