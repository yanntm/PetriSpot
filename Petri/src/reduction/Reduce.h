#pragma once
#include "reduction/Coordinator.h"

namespace petri::reduction {
template<class T> struct Result {
  SparsePetriNet<T> net;
  std::vector<size_t> placeMap, transitionMap;
  std::array<RuleStats, 11> stats;
  bool limited = false;
  size_t edits = 0;
};

/** Input has no external metadata contract. STATESPACE conservatively retains
 * token coordinates and duplicate transitions. Trace/image reconstruction and
 * arbitrary external records are not accepted by this interface yet. */
template<class T, class Trace>
Result<T> reduce(SparsePetriNet<T> net, Configuration config,
                 std::vector<bool> observed, Trace& trace) {
  Workspace<T> workspace(std::move(net), config, std::move(observed));
  Coordinator<T, Trace> coordinator(workspace, trace);
  coordinator.execute();
  Result<T> result;
  result.net = workspace.publish(result.placeMap, result.transitionMap);
  result.stats = coordinator.stats; result.limited = workspace.limited; result.edits = workspace.changes;
  return result;
}

template<class T>
Result<T> reduce(SparsePetriNet<T> net, Configuration config, std::vector<bool> observed = {}) {
  NoTrace trace;
  return reduce(std::move(net), config, std::move(observed), trace);
}
}
