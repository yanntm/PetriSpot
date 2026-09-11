#pragma once
#include <ostream>
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
 * otherwise (Counting.h). Any other goal ignores it. `knownDead` are
 * transitions a caller proved never enabled (an outside test): retired
 * first, as dead, so the rules start from what that proof exposes.
 * Trace/image reconstruction is not accepted by this interface yet. */
template<class T, class Trace>
Result<T> reduce(SparsePetriNet<T> net, Configuration config,
                 std::vector<bool> observed, Trace& trace,
                 std::optional<Counting<T>> counting = std::nullopt,
                 std::vector<size_t> knownDead = {}) {
  Workspace<T> workspace(std::move(net), config, std::move(observed));
  for (size_t t : knownDead) if (t >= workspace.transitions.size()) throw std::invalid_argument("Known dead transition out of range");
  workspace.retireDeadTransitions(std::move(knownDead));
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

/** The dead-transition statistics of a result as one diagnostics line. */
inline void describeDeadTransitions(const DeadStats& d, std::ostream& os) {
  if (d.tested == 0 && d.places == 0 && !d.limited) return;
  os << "Reduction dead transitions: " << d.found << " (" << d.byBound << " by the bounds of " << d.places
      << " places, " << d.stuck << " never marked, " << d.placesMs << " ms), then " << d.tested << " tested, " << d.solves << " solves, "
      << d.pivots << " pivots, " << d.passes << " passes, " << d.spentMs << " ms" << (d.limited ? " (budget)" : "") << ".\n";
}

template<class T>
Result<T> reduce(SparsePetriNet<T> net, Configuration config, std::vector<bool> observed = {},
                 std::optional<Counting<T>> counting = std::nullopt, std::vector<size_t> knownDead = {}) {
  NoTrace trace;
  return reduce(std::move(net), config, std::move(observed), trace, std::move(counting), std::move(knownDead));
}
}
