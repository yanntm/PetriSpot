#pragma once
#include <optional>
#include <ostream>
#include "expr/Property.h"
#include "reduction/Reduce.h"

namespace petri::reduction {
inline void collectSupport(const expr::Expression& expression, std::vector<bool>& support) {
  for (const auto& [p, coefficient] : expression.atom.terms) {
    (void) coefficient;
    if (p >= support.size()) throw std::invalid_argument("Property place outside reduction input");
    support[p] = true;
  }
  for (const auto& child : expression.children) collectSupport(child, support);
}
inline void collectSupport(const expr::CtlFormula& formula, std::vector<bool>& support) {
  collectSupport(formula.pred, support);
  for (const auto& child : formula.kids) collectSupport(child, support);
}
inline void remap(expr::Expression& expression, const std::vector<size_t>& map) {
  for (auto& term : expression.atom.terms) {
    if (term.first >= map.size() || map[term.first] == std::numeric_limits<size_t>::max())
      throw std::logic_error("Reduction lost an observed place");
    term.first = map[term.first];
  }
  for (auto& child : expression.children) remap(child, map);
}
inline void remap(expr::CtlFormula& formula, const std::vector<size_t>& map) {
  remap(formula.pred, map);
  for (auto& child : formula.kids) remap(child, map);
}

/** Derive a joint goal. General CTL and mixed deadlock/state queries use the
 * strong local subset (LTL configuration); no stutter assumption is inferred.
 * Hints and executable witnesses retain the original net until lifting exists.
 * This adapter owns diagnostics; the reduction kernel writes no output. */
template<class T>
std::optional<Result<T>> prepareQueries(const SparsePetriNet<T>& original,
    std::vector<expr::Property>& properties, bool enabled, bool needsOriginal,
    Configuration config, std::ostream& diagnostics) {
  if (!enabled) return std::nullopt;
  for (const auto& p : properties) needsOriginal |= !p.hint.empty();
  if (needsOriginal) {
    diagnostics << "Reduction skipped: original trace/hint coordinates requested.\n";
    return std::nullopt;
  }
  bool temporal = false, deadlock = false, state = false;
  std::vector<bool> support(original.getPlaceCount(), false);
  for (const auto& p : properties) {
    if (p.kind == expr::PropertyKind::Unsupported) {
      diagnostics << "Reduction skipped: unsupported property in query set.\n";
      return std::nullopt;
    }
    temporal |= p.kind == expr::PropertyKind::CTL;
    deadlock |= p.kind == expr::PropertyKind::Deadlock;
    state |= p.kind != expr::PropertyKind::CTL && p.kind != expr::PropertyKind::Deadlock;
    collectSupport(p.body, support); collectSupport(p.ctl, support);
  }
  config.goal = temporal || (deadlock && state) ? Goal::LTL
      : deadlock ? Goal::DEADLOCK : Goal::REACHABILITY;
  auto start = std::chrono::steady_clock::now();
  auto result = reduce(original, config, std::move(support));
  for (auto& property : properties) { remap(property.body, result.placeMap); remap(property.ctl, result.placeMap); }
  diagnostics << "Reduction (native subset): " << original.getPlaceCount() << " -> "
      << result.net.getPlaceCount() << " places, " << original.getTransitionCount() << " -> "
      << result.net.getTransitionCount() << " transitions, " << original.getArcCount() << " -> "
      << result.net.getArcCount() << " arcs, "
      << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
      << " ms" << (result.limited ? " (limit)" : "") << ".\n";
  for (const auto& stat : result.stats) if (stat.edits)
    diagnostics << "Reduction rule " << stat.name << ": " << stat.edits << " sparse edits.\n";
  return result;
}
}
