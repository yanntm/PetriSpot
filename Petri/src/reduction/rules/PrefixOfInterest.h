#pragma once
#include "reduction/graph/Dependency.h"
#include "reduction/graph/Stabilizing.h"
namespace petri::reduction {
/** StructuralReduction.findSCCSuffixes/computeSafeNodes. Reachability starts at
 * property support; temporal modes also keep cyclic behavior. Deadlock starts
 * only at cyclic SCCs after excluding stabilizing transitions. All modes then
 * keep predecessors; temporal modes add shared consumption guards and reclose.
 * Deadlock deductions are results, never encoded by destructive net edits. */
struct PrefixOfInterest {
  static constexpr const char* name = "prefix-of-interest";
  template<class T> static void apply(Workspace<T>& w) {
    Goal goal = w.config.goal;
    if (!w.config.relevance || goal == Goal::LTL || goal == Goal::LIVENESS || goal == Goal::STATESPACE) return;
    bool deadlock = goal == Goal::DEADLOCK, temporal = temporalPrefix(goal);
    if (deadlock) {
      // Java also proves this in ruleReduceTrans. Do not let a graph that has
      // no input edges for sources hide an unconditional enabled transition.
      for (size_t t = 0; t < w.transitions.size(); ++t)
        if (w.liveT[t] && w.pre.getColumn(t).size() == 0) { w.deadlock = false; return; }
    }
    std::vector<bool> skipped;
    if (deadlock) skipped = computeStabilizing(w).transitions;
    if (w.stop()) return;
    Graph graph = dependencyGraph(w, skipped);
    if (w.stop()) return;
    std::vector<bool> kept(w.places.size());
    if (deadlock || temporal) {
      auto components = graph.components([&] { return w.stop(); });
      if (w.stop()) return;
      bool cyclic = false;
      for (const auto& component : components) {
        size_t p = component.front();
        if (component.size() == 1 && !std::binary_search(graph.next[p].begin(), graph.next[p].end(), p)) continue;
        cyclic = true;
        for (size_t q : component) kept[q] = true;
      }
      if (deadlock && !cyclic) { w.deadlock = true; return; }
    }
    if (!deadlock) observationSeeds(w, kept);
    if (deadlock) graph = dependencyGraph(w);
    graph.prefix(kept, [&] { return w.stop(); });
    if (w.stop()) return;
    if (temporal) {
      // Deliberately preserve Java's transition-order scan with immediate seed
      // updates, followed by one complete prefix closure.
      for (size_t t = 0; t < w.transitions.size(); ++t) if (w.liveT[t]) {
        const auto& pre = w.pre.getColumn(t);
        bool touches = false;
        for (size_t i = 0; i < pre.size(); ++i) if (kept[pre.keyAt(i)]) { touches = true; break; }
        if (touches) for (size_t i = 0; i < pre.size(); ++i) kept[pre.keyAt(i)] = true;
      }
      graph.prefix(kept, [&] { return w.stop(); });
    }
    if (!w.stop()) dropOutsidePrefix(w, kept);
  }
};
}
