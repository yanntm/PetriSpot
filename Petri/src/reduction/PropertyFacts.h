#pragma once
#include "core/Arithmetic.hpp"
#include "core/SparsePetriNet.h"
#include "expr/CtlSimplify.h"
#include "expr/Property.h"

namespace petri::reduction {

/** Facts available without exploration. A place is constant when no transition
 * changes it, or when it lies in the greatest initially empty siphon (it never
 * receives its first token). Initial deadlock makes every place constant. A
 * transition enabled using only constant guards rules out every deadlock. */
struct PropertyFacts {
  std::vector<bool> constant;
  bool initialDeadlock = true;
  bool noDeadlocks = false;
};

/** The greatest siphon of initially empty places: start from them all, drop the
 * outputs of every transition none of whose inputs is a candidate, propagate. */
template<class T> std::vector<bool> emptySiphon(const SparsePetriNet<T>& net) {
  size_t places = net.getPlaceCount(), transitions = net.getTransitionCount();
  std::vector<bool> candidate(places);
  for (size_t p = 0; p < places; ++p) candidate[p] = net.getMarks()[p] == 0;
  std::vector<size_t> inputs(transitions, 0), todo;
  for (size_t t = 0; t < transitions; ++t) {
    const auto& pre = net.getFlowPT().getColumn(t);
    for (size_t i = 0; i < pre.size(); ++i) inputs[t] += candidate[pre.keyAt(i)] ? 1 : 0;
    if (inputs[t] == 0) todo.push_back(t);
  }
  auto consumers = net.getFlowPT().transpose();
  for (size_t k = 0; k < todo.size(); ++k) {
    const auto& post = net.getFlowTP().getColumn(todo[k]);
    for (size_t i = 0; i < post.size(); ++i) {
      size_t p = post.keyAt(i);
      if (!candidate[p]) continue;
      candidate[p] = false;
      const auto& row = consumers.getColumn(p);
      for (size_t j = 0; j < row.size(); ++j) if (--inputs[row.keyAt(j)] == 0) todo.push_back(row.keyAt(j));
    }
  }
  return candidate;
}

template<class T> PropertyFacts propertyFacts(const SparsePetriNet<T>& net) {
  PropertyFacts facts;
  facts.constant.assign(net.getPlaceCount(), true);
  for (size_t t = 0; t < net.getTransitionCount(); ++t) {
    const auto& pre = net.getFlowPT().getColumn(t);
    const auto& post = net.getFlowTP().getColumn(t);
    bool enabled = true;
    for (size_t i = 0; i < pre.size(); ++i)
      if (net.getMarks()[pre.keyAt(i)] < pre.valueAt(i)) enabled = false;
    facts.initialDeadlock &= !enabled;
    for (const auto* col : {&pre, &post})
      for (size_t i = 0; i < col->size(); ++i) {
        size_t p = col->keyAt(i);
        if (pre.get(p) != post.get(p)) facts.constant[p] = false;
      }
  }
  if (facts.initialDeadlock) {
    std::fill(facts.constant.begin(), facts.constant.end(), true);
    return facts;
  }
  auto siphon = emptySiphon(net);
  for (size_t p = 0; p < siphon.size(); ++p) if (siphon[p]) facts.constant[p] = true;
  for (size_t t = 0; t < net.getTransitionCount(); ++t) {
    const auto& pre = net.getFlowPT().getColumn(t);
    bool always = true;
    for (size_t i = 0; i < pre.size(); ++i) {
      size_t p = pre.keyAt(i);
      if (!facts.constant[p] || net.getMarks()[p] < pre.valueAt(i)) { always = false; break; }
    }
    if (always) { facts.noDeadlocks = true; break; }
  }
  return facts;
}

inline long long subtractConstant(long long value, long long offset) {
  if ((offset > 0 && value < std::numeric_limits<long long>::min() + offset)
      || (offset < 0 && value > std::numeric_limits<long long>::max() + offset))
    throw std::overflow_error("Property constant substitution overflow");
  return value - offset;
}

/** Remove constant terms from an atomic proposition before boolean
 * simplification. If its existing long-long representation cannot express the
 * rewritten atom, leave that atom intact rather than narrowing a marking. */
template<class T>
expr::Expression substituteConstants(expr::Expression expression, const PropertyFacts& facts,
                                     const std::vector<T>& marking) {
  if (expression.kind == expr::Expression::Kind::Atom) {
    auto atom = expression.atom;
    atom.terms.clear();
    long long offset = 0;
    try {
      for (const auto& [p, coefficient] : expression.atom.terms) {
        if (facts.constant.at(p))
          offset = petri::addExact(offset, petri::multiplyExact(coefficient, petri::castExact<long long>(marking.at(p))));
        else atom.terms.emplace_back(p, coefficient);
      }
      atom.constant = subtractConstant(atom.constant, offset);
    } catch (const std::overflow_error&) { return expression; }
    expression.atom = std::move(atom);
  }
  for (auto& child : expression.children) child = substituteConstants(std::move(child), facts, marking);
  return expr::simplify(expression);
}

template<class T>
expr::CtlFormula substituteConstants(expr::CtlFormula formula, const PropertyFacts& facts,
                                    const std::vector<T>& marking) {
  using expr::CtlOp;
  formula.pred = substituteConstants(std::move(formula.pred), facts, marking);
  for (auto& child : formula.kids) child = substituteConstants(std::move(child), facts, marking);
  formula = expr::ctlSimplify(formula);
  // Simplification can introduce a deadlock leaf, e.g. EX true. Substitute
  // that fact as well, so the enclosing expression sees a boolean constant.
  if ((formula.op == CtlOp::Deadlock || formula.op == CtlOp::NoDeadlock)
      && (facts.initialDeadlock || facts.noDeadlocks))
    return expr::CtlFormula::constant((formula.op == CtlOp::Deadlock) == facts.initialDeadlock);
  return formula;
}

/** Propagate structural facts into the normal property representation. Boolean
 * formulas become constants through the existing simplifiers. Constant bound
 * forms retain their coordinates (the syntax has no affine offset), with an
 * exact known bound so normal initial-goal evaluation can close them. */
template<class T>
size_t simplifyProperties(const SparsePetriNet<T>& net, std::vector<expr::Property>& properties) {
  const auto facts = propertyFacts(net);
  size_t changed = 0;
  for (auto& property : properties) {
    const expr::Property before = property;
    if (property.kind == expr::PropertyKind::Bound) {
      bool constant = true;
      long long value = 0;
      try {
        for (const auto& [p, coefficient] : property.boundForm().terms) {
          if (!facts.constant.at(p)) { constant = false; break; }
          value = petri::addExact(value, petri::multiplyExact(coefficient, petri::castExact<long long>(net.getMarks().at(p))));
        }
        if (constant && value >= 0) property.boundHint = value;
      } catch (const std::overflow_error&) { /* Keep the existing representable form. */ }
    } else if (property.kind == expr::PropertyKind::Deadlock) {
      if (facts.initialDeadlock || facts.noDeadlocks) {
        property.kind = expr::PropertyKind::Reachability;
        property.body = expr::Expression::constant(facts.initialDeadlock);
      }
    } else if (property.kind == expr::PropertyKind::CTL) {
      property.ctl = substituteConstants(std::move(property.ctl), facts, net.getMarks());
    } else if (property.kind != expr::PropertyKind::Unsupported) {
      property.body = substituteConstants(std::move(property.body), facts, net.getMarks());
    }
    if (property.kind != before.kind || !(property.body == before.body) || !(property.ctl == before.ctl)
        || property.boundHint != before.boundHint) ++changed;
  }
  return changed;
}
}
