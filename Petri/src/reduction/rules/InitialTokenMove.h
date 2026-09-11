#pragma once
#include "core/Arithmetic.hpp"
#include "reduction/Workspace.h"
namespace petri::reduction {
/** Java ruleReducePlaces(moveTokens): at final stability, a marked unobserved
 * source place controls one invisible transition by itself. Fire its complete
 * continuation and discard any unusable remainder. This changes the initial
 * state only through unobserved steps; reachability/deadlock are preserved. */
struct InitialTokenMove {
  static constexpr const char* name = "initial-token-move";
  template<class T> static void apply(Workspace<T>& w) {
    if (w.config.goal != Goal::REACHABILITY && w.config.goal != Goal::DEADLOCK) return;
    for (size_t p = w.places.size(); p-- > 0;) {
      if (w.stop()) return;
      if (!w.liveP[p] || w.observed[p] || w.marks[p] == 0 || w.producers.getColumn(p).size() != 0) continue;
      const auto& row = w.consumers.getColumn(p);
      if (row.size() != 1) continue;
      size_t t = row.keyAt(0);
      if (w.pre.getColumn(t).size() != 1 || w.visible(t)) continue;
      T count = w.marks[p] / row.valueAt(0);
      SparseArray<T> updated;
      const auto& post = w.post.getColumn(t);
      for (size_t i = 0; i < post.size(); ++i) {
        size_t q = post.keyAt(i);
        updated.put(q, petri::addExact(w.marks[q], petri::multiplyExact(count, post.valueAt(i))));
      }
      w.marks[p] = 0;
      for (size_t i = 0; i < updated.size(); ++i) w.marks[updated.keyAt(i)] = updated.valueAt(i);
      ++w.changes;
    }
  }
};
}
