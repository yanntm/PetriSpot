/*
 * NetStats.h
 *
 * Structural statistics of a WalkNet, printed as histograms: transition
 * arities (pre, post, effect sizes), place fan-out (consumers), tokens.
 */
#ifndef PETRI_WALK_NETSTATS_H_
#define PETRI_WALK_NETSTATS_H_

#include <map>
#include <ostream>

#include "walk/WalkNet.h"

namespace petri::walk
{

namespace detail
{
inline void printHistogram (std::ostream &os, const char *title,
                            const std::map<size_t, size_t> &h)
{
  os << title << " :";
  for (const auto &kv : h) os << " " << kv.first << "->" << kv.second;
  os << std::endl;
}
} // namespace detail

template<typename T>
  void printNetStats (std::ostream &os, const WalkNet<T> &net)
  {
    std::map<size_t, size_t> pre, post, effect, fanout;
    size_t emptyEffect = 0;
    for (size_t t = 0; t < net.transitionCount (); ++t) {
      ++pre[net.pre (t).size ()];
      ++post[net.getNet ().getFlowTP ().getColumn (t).size ()];
      ++effect[net.effect (t).size ()];
      if (!net.hasEffect (t)) ++emptyEffect;
    }
    size_t tokens = 0, marked = 0;
    for (size_t p = 0; p < net.placeCount (); ++p) {
      ++fanout[net.consumersOf (p).size ()];
      T m = net.initialMarking ()[p];
      if (m != 0) {
        ++marked;
        tokens += static_cast<size_t> (m);
      }
    }
    os << "Net: " << net.placeCount () << " places, " << net.transitionCount ()
        << " transitions, " << emptyEffect << " with empty effect; initial marking has "
        << tokens << " tokens on " << marked << " places." << std::endl;
    detail::printHistogram (os, "pre-size -> #transitions", pre);
    detail::printHistogram (os, "post-size -> #transitions", post);
    detail::printHistogram (os, "effect-size -> #transitions", effect);
    detail::printHistogram (os, "fan-out -> #places", fanout);
  }

} // namespace petri::walk

#endif /* PETRI_WALK_NETSTATS_H_ */
