#pragma once
#include "reduction/Coordinator.h"

namespace petri::reduction {
template<class T> struct Result {
  SparsePetriNet<T> net;
  std::vector<size_t> placeMap, transitionMap;
  std::array<RuleStats, 27> stats;
  bool limited = false;
  std::optional<bool> deadlock;
  std::optional<Counting<T>> counting; // the published net's record, STATESPACE only
  DeadStats dead;
  size_t edits = 0;
};

/** `counting` is the record of the input net, attached under STATESPACE only:
 * the identity when the input vouches for itself, the blocks it carried
 * otherwise (Counting.h). Any other goal ignores it. Trace/image
 * reconstruction is not accepted by this interface yet. */
template<class T, class Trace>
Result<T> reduce(SparsePetriNet<T> net, Configuration config,
                 std::vector<bool> observed, Trace& trace,
                 std::optional<Counting<T>> counting = std::nullopt) {
  Workspace<T> workspace(std::move(net), config, std::move(observed));
  if (config.goal == Goal::STATESPACE) {
    if (!counting) counting = Counting<T>::identity(workspace.places.size(), workspace.transitions.size());
    if (counting->pcoef.size() != workspace.places.size()
        || (counting->tmult && counting->tmult->size() != workspace.transitions.size()))
      throw std::invalid_argument("Counting record does not fit the net");
    workspace.counting = std::move(counting);
  }
  Coordinator<T, Trace> coordinator(workspace, trace);
  coordinator.execute();
  Result<T> result;
  result.net = workspace.publish(result.placeMap, result.transitionMap);
  if (workspace.counting) result.counting = workspace.counting->compact(result.placeMap, result.transitionMap);
  result.deadlock = workspace.deadlock; result.dead = workspace.dead;
  result.stats = coordinator.stats; result.limited = workspace.limited; result.edits = workspace.changes;
  return result;
}

template<class T>
Result<T> reduce(SparsePetriNet<T> net, Configuration config, std::vector<bool> observed = {},
                 std::optional<Counting<T>> counting = std::nullopt) {
  NoTrace trace;
  return reduce(std::move(net), config, std::move(observed), trace, std::move(counting));
}
}
