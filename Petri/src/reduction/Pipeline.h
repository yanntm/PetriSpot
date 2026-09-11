#pragma once
#include <optional>
#include <ostream>
#include "expr/InitialState.h"
#include "reduction/Properties.h"
#include "reduction/cli/PropertyResults.h"

namespace petri::reduction {
/** The net the engines analyse and what produced it. `net` is the original
 * when no reduction ran; `reduction` is the last reduction result otherwise. */
template<class T> struct Prepared {
  SparsePetriNet<T> net;
  std::optional<Result<T>> reduction;
  size_t rounds = 0;
};

/** The cheap preparation every engine runs before spending anything, to a
 * fixpoint: the constant places of the current net are substituted into the
 * formulas and the formulas simplified; what the initial marking decides is
 * answered and an until it decides rewritten, a formula fallen into the
 * reachability fragment requalified; decided properties are reported and
 * dropped, so the support the reduction protects shrinks; then the net is
 * reduced for the remaining kinds and supports. A reduction that edited the
 * net or a property that changed starts another round, since each opens the
 * other's way; the loop ends at stability, when no property is left, or at
 * `config.maxRounds`. With `reduce` off the net stays as it is and only the
 * formula side runs: the degraded mode, still worth having.
 *
 * `answers` receives the FORMULA lines of the decided properties; null keeps
 * them in `properties` as constants instead (the standalone export). Before
 * any reduction the technique is INITIAL_STATE, after one STRUCTURAL_REDUCTION. */
template<class T>
Prepared<T> prepare(const SparsePetriNet<T>& original, std::vector<expr::Property>& properties,
                    bool reduce, bool needsOriginal, Configuration config,
                    std::ostream* answers, std::ostream& diagnostics) {
  Prepared<T> out;
  out.net = original;
  const char* technique = "TOPOLOGICAL INITIAL_STATE";
  for (size_t round = 0;; ++round) {
    bool changed = simplifyProperties(out.net, properties) > 0;
    changed |= expr::simplifyInitial(out.net, properties) > 0;
    if (answers) consumeSolvedProperties(out.net, properties, *answers, technique);
    if (properties.empty() || !reduce || round >= config.maxRounds) break;
    auto result = prepareQueries(out.net, properties, true, needsOriginal, config, diagnostics);
    if (!result) break;
    bool edited = result->edits > 0 || result->deadlock.has_value();
    if (edited) {
      out.net = std::move(result->net);
      result->net = SparsePetriNet<T>();
      out.reduction = std::move(*result);
      technique = "TOPOLOGICAL STRUCTURAL_REDUCTION";
      ++out.rounds;
    }
    if (!edited && !changed) break;
  }
  return out;
}
}
