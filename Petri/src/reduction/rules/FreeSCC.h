#pragma once
#include "reduction/graph/Graph.h"
#include "core/Arithmetic.hpp"
#include "reduction/Workspace.h"
namespace petri::reduction {
/** StructuralReduction.findFreeSCC: unobserved places joined by unit one-input,
 * one-output transitions can freely redistribute their total tokens. Fuse each
 * nontrivial SCC by summing markings and arc weights. Internal transfers become
 * self loops, preserving their ability to stutter or prevent a deadlock.
 * LTL forbids fusion; STATESPACE needs the reference's counting records. */
struct FreeSCC {
  static constexpr const char* name = "free-scc";
  template<class T> static void apply(Workspace<T>& w) {
    if (w.config.goal == Goal::LTL || w.config.goal == Goal::STATESPACE) return;
    Graph graph(w.places.size());
    for (size_t t = 0; t < w.transitions.size(); ++t) if (w.liveT[t]) {
      if (w.stop()) return;
      const auto& pre = w.pre.getColumn(t); const auto& post = w.post.getColumn(t);
      if (pre.size() == 1 && post.size() == 1 && pre.valueAt(0) == 1 && post.valueAt(0) == 1
          && !w.observed[pre.keyAt(0)] && !w.observed[post.keyAt(0)]) graph.edge(pre.keyAt(0), post.keyAt(0));
    }
    graph.finish();
    auto components = graph.components([&] { return w.stop(); });
    if (w.stop()) return;
    for (const auto& component : components) if (component.size() > 1) {
      size_t kept = component.front();
      T marking = w.marks[kept];
      SparseArray<T> inputs = w.consumers.getColumn(kept), outputs = w.producers.getColumn(kept);
      // Check all sums before modifying this component.
      for (size_t i = 1; i < component.size(); ++i) {
        size_t p = component[i];
        marking = petri::addExact(marking, w.marks[p]);
        for (auto pair : {std::pair{&inputs, &w.consumers.getColumn(p)}, std::pair{&outputs, &w.producers.getColumn(p)}})
          for (size_t j = 0; j < pair.second->size(); ++j) {
            size_t t = pair.second->keyAt(j);
            pair.first->put(t, petri::addExact(pair.first->get(t), pair.second->valueAt(j)));
          }
      }
      if (w.stop()) return;
      w.safe = false; // the fused place holds the component's total
      for (size_t i = 1; i < component.size(); ++i) w.retirePlace(component[i]);
      w.marks[kept] = marking;
      for (size_t i = 0; i < inputs.size(); ++i) {
        size_t t = inputs.keyAt(i); auto col = w.pre.getColumn(t); col.put(kept, inputs.valueAt(i));
        w.replacePre(t, std::move(col));
      }
      for (size_t i = 0; i < outputs.size(); ++i) {
        size_t t = outputs.keyAt(i); auto col = w.post.getColumn(t); col.put(kept, outputs.valueAt(i));
        w.replacePost(t, std::move(col));
      }
    }
  }
};
}
