#pragma once
#include <algorithm>
#include <optional>
#include <ostream>
#include "reduction/PropertyFacts.h"

namespace petri::reduction {
/** Read a result from a simplified body, without running an analysis. Bounds
 * are closed when the initial value attains their known upper bound. */
template<class T>
std::optional<std::string> propertyResult(const SparsePetriNet<T>& net, const expr::Property& property) {
  using expr::PropertyKind;
  if ((property.kind == PropertyKind::Reachability || property.kind == PropertyKind::Invariant)
      && property.body.isConstant())
    return property.body.kind == expr::Expression::Kind::True ? "TRUE" : "FALSE";
  if (property.kind == PropertyKind::CTL && property.ctl.isConstant())
    return property.ctl.isTrue() ? "TRUE" : "FALSE";
  if (property.kind == PropertyKind::Bound && property.boundHint >= 0) {
    long long initial = 0;
    try {
      for (const auto& [p, coefficient] : property.boundForm().terms)
        initial = petri::addExact(initial, petri::multiplyExact(coefficient,
            petri::castExact<long long>(net.getMarks().at(p))));
      if (initial == property.boundHint) return std::to_string(initial);
    } catch (const std::overflow_error&) { /* Leave evaluation to the consumer. */ }
  }
  return std::nullopt;
}

/** Consume decided properties before allocating walk, CTL or LP machinery. */
template<class T>
size_t consumeSolvedProperties(const SparsePetriNet<T>& net, std::vector<expr::Property>& properties,
                               std::ostream& output) {
  return std::erase_if(properties, [&](const auto& property) {
    auto result = propertyResult(net, property);
    if (!result) return false;
    output << "FORMULA " << property.name << " " << *result
        << " TECHNIQUES TOPOLOGICAL STRUCTURAL_REDUCTION\n";
    if (property.kind == expr::PropertyKind::Bound)
      output << "BOUND " << property.name << " " << *result << '\n';
    return true;
  });
}
}
