#pragma once
#include <unordered_map>
#include "reduction/Workspace.h"

namespace petri::reduction {
/** Identical incidence rows have a fixed marking difference. The smaller
 * initial marking entails all guards of the larger one. Keep both when the
 * larger coordinate is observed. Never remove a previously retained guard
 * without a live representative. STATESPACE retains token-count coordinates. */
struct DuplicatePlace {
  static constexpr const char* name = "duplicate-place";
  template<class T> static void apply(Workspace<T>& w) {
    if (w.config.goal == Goal::STATESPACE) return;
    std::unordered_map<size_t, std::vector<size_t>> buckets;
    for (size_t p = 0; p < w.places.size(); ++p) if (w.liveP[p]) {
      if (p % 256 == 0 && w.stop()) return;
      const auto& in = w.producers.getColumn(p); const auto& out = w.consumers.getColumn(p);
      auto& bucket = buckets[in.hash() * 31 + out.hash()];
      bool matched = false;
      for (size_t& other : bucket) {
        if (!(in == w.producers.getColumn(other)) || !(out == w.consumers.getColumn(other))) continue;
        if (w.marks[p] >= w.marks[other]) {
          if (!w.observed[p]) w.retirePlace(p);
        } else {
          if (!w.observed[other]) w.retirePlace(other);
          other = p;
        }
        matched = true; break;
      }
      if (!matched) bucket.push_back(p);
    }
  }
};
}
