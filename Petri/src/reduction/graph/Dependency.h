#pragma once
#include "reduction/graph/Graph.h"
#include "reduction/Workspace.h"
namespace petri::reduction {
inline bool temporalPrefix(Goal goal) {
  return goal == Goal::SI_LTL || goal == Goal::LI_LTL || goal == Goal::SI_CTL;
}
/** Java buildGraph represented in the forward direction; previous is the
 * adjacency used by collectPrefix. The skipped mask affects edges, not support. */
template<class T> Graph dependencyGraph(Workspace<T>& w, const std::vector<bool>& skipped = {}) {
  Graph graph(w.places.size());
  bool all = temporalPrefix(w.config.goal) || w.config.goal == Goal::DEADLOCK;
  for (size_t t = 0; t < w.transitions.size(); ++t) if (w.liveT[t] && (skipped.empty() || !skipped[t])) {
    if (w.stop()) return graph;
    const auto& pre = w.pre.getColumn(t); const auto& post = w.post.getColumn(t);
    for (size_t j = 0; j < post.size(); ++j) {
      size_t q = post.keyAt(j);
      if (!all && pre.get(q) == post.valueAt(j)) continue;
      for (size_t i = 0; i < pre.size(); ++i) if (all || pre.keyAt(i) != q) graph.edge(pre.keyAt(i), q);
    }
  }
  graph.finish(); return graph;
}
/** Property seeds from computeSafeNodes, including consuming-only visible
 * effects. This is separate from the graph's output-edge predicate. */
template<class T> void observationSeeds(const Workspace<T>& w, std::vector<bool>& kept) {
  for (size_t p = 0; p < kept.size(); ++p) kept[p] = kept[p] || w.observed[p];
  for (size_t t = 0; t < w.transitions.size(); ++t) if (w.liveT[t] && w.visible(t)) {
    const auto& pre = w.pre.getColumn(t);
    for (size_t i = 0; i < pre.size(); ++i) kept[pre.keyAt(i)] = true;
  }
}
/** Exact dropIrrelevant/dropPlaces(andOutputs) behavior: consumers losing their
 * entire preset disappear; consumers retaining a guard remain. */
template<class T> void dropOutsidePrefix(Workspace<T>& w, const std::vector<bool>& kept) {
  std::vector<bool> consumers(w.transitions.size());
  for (size_t p = 0; p < kept.size(); ++p) if (w.liveP[p] && !kept[p]) {
    const auto& row = w.consumers.getColumn(p);
    for (size_t i = 0; i < row.size(); ++i) consumers[row.keyAt(i)] = true;
    w.retirePlace(p);
  }
  for (size_t t = 0; t < consumers.size(); ++t)
    if (consumers[t] && w.liveT[t] && w.pre.getColumn(t).size() == 0) w.retireTransition(t);
}
}
