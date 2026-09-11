#pragma once
#include "reduction/rules/ConstantPlace.h"

namespace petri::reduction {
/** Greatest initially empty siphon. Remove candidate outputs of every
 * transition with no candidate input; propagate lost inputs with counters.
 * Remaining places cannot acquire their first token, so their consumers are
 * dead. O(P+T+arcs), including source transitions. LIVENESS skips removal;
 * STATESPACE retains zero place coordinates. */
struct EmptySiphon {
  static constexpr const char* name = "empty-siphon";
  template<class T> static void apply(Workspace<T>& w) {
    if (w.config.goal == Goal::LIVENESS) return;
    std::vector<bool> candidate(w.places.size());
    std::vector<size_t> inputs(w.transitions.size(), 0), todo;
    for (size_t p = 0; p < candidate.size(); ++p) candidate[p] = w.liveP[p] && w.marks[p] == 0;
    for (size_t t = 0; t < inputs.size(); ++t) if (w.liveT[t]) {
      if (t % 256 == 0 && w.stop()) return;
      const auto& pre = w.pre.getColumn(t);
      for (size_t i = 0; i < pre.size(); ++i) inputs[t] += candidate[pre.keyAt(i)] ? 1 : 0;
      if (inputs[t] == 0) todo.push_back(t);
    }
    for (size_t k = 0; k < todo.size(); ++k) {
      if (k % 256 == 0 && w.stop()) return;
      const auto& post = w.post.getColumn(todo[k]);
      for (size_t i = 0; i < post.size(); ++i) {
        size_t p = post.keyAt(i);
        if (!candidate[p]) continue;
        candidate[p] = false;
        const auto& consumers = w.consumers.getColumn(p);
        for (size_t j = 0; j < consumers.size(); ++j) {
          size_t t = consumers.keyAt(j);
          if (--inputs[t] == 0) todo.push_back(t);
        }
      }
    }
    for (size_t p = 0; p < candidate.size(); ++p) if (candidate[p]) {
      const auto consumers = w.consumers.getColumn(p);
      for (size_t i = 0; i < consumers.size(); ++i) w.retireTransition(consumers.keyAt(i));
      if (w.observed[p] || w.config.goal == Goal::STATESPACE) w.erasePlaceArcs(p); else w.retirePlace(p);
    }
  }
};
}
