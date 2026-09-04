/*
 * NetResolver.h
 *
 * Name-to-index lookup over a net and the desugaring of fireability into
 * marking conditions, shared by the property parsers.
 */
#ifndef PETRI_PARSE_NETRESOLVER_H_
#define PETRI_PARSE_NETRESOLVER_H_

#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "core/SparsePetriNet.h"
#include "expr/Expression.h"

namespace petri
{

template<typename T>
  class NetResolver
  {
    const SparsePetriNet<T> &net;
    std::unordered_map<std::string, size_t> placeIndex;
    std::unordered_map<std::string, size_t> transitionIndex;

  public:
    explicit NetResolver (const SparsePetriNet<T> &n)
        : net (n)
    {
      const auto &pn = net.getPnames ();
      for (size_t i = 0; i < pn.size (); ++i) placeIndex[pn[i]] = i;
      const auto &tn = net.getTnames ();
      for (size_t i = 0; i < tn.size (); ++i) transitionIndex[tn[i]] = i;
    }

    const SparsePetriNet<T>& getNet () const
    {
      return net;
    }
    size_t placeCount () const
    {
      return net.getPlaceCount ();
    }
    size_t transitionCount () const
    {
      return net.getTransitionCount ();
    }
    std::optional<size_t> place (const std::string &name) const
    {
      auto it = placeIndex.find (name);
      if (it == placeIndex.end ()) return std::nullopt;
      return it->second;
    }
    std::optional<size_t> transition (const std::string &name) const
    {
      auto it = transitionIndex.find (name);
      if (it == transitionIndex.end ()) return std::nullopt;
      return it->second;
    }

    /** is-fireable(t): every pre-place holds at least the arc weight. */
    petri::expr::Expression fireable (size_t t) const
    {
      using petri::expr::Expression;
      const auto &pre = net.getFlowPT ().getColumn (t);
      std::vector<Expression> kids;
      kids.reserve (pre.size ());
      for (size_t i = 0; i < pre.size (); ++i) {
        petri::expr::LinearAtom a;
        a.addTerm (pre.keyAt (i), 1);
        a.constant = static_cast<long long> (pre.valueAt (i));
        a.op = petri::expr::Cmp::GE;
        kids.push_back (Expression::makeAtom (std::move (a)));
      }
      if (kids.empty ()) return Expression::constant (true);
      if (kids.size () == 1) return std::move (kids[0]);
      return Expression::makeAnd (std::move (kids));
    }

    /** Some transition of the list is fireable; false for an empty list. */
    petri::expr::Expression anyFireable (const std::vector<size_t> &ts) const
    {
      using petri::expr::Expression;
      std::vector<Expression> kids;
      kids.reserve (ts.size ());
      for (size_t t : ts) kids.push_back (fireable (t));
      if (kids.empty ()) return Expression::constant (false);
      if (kids.size () == 1) return std::move (kids[0]);
      return Expression::makeOr (std::move (kids));
    }
  };

} // namespace petri

#endif /* PETRI_PARSE_NETRESOLVER_H_ */
