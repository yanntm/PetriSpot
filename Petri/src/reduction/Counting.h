#pragma once
#include <limits>
#include <optional>
#include <string>
#include <vector>
#include "core/Arithmetic.hpp"

namespace petri::reduction {
/** What each object of the net being edited stands for in the net it came
 * from, kept only when someone will count. Weights are stored minus one, so
 * the identity is all zeroes. `tmult` is absent once the arcs of this net are
 * no longer those of the baseline; `arcsLost` then says why. `pcoef` and
 * `pdrop` stay valid through every rule that may run while the record is
 * attached (see algorithm.md, section 4). Indexed by workspace slot until
 * `compact` renumbers to a published net. */
template<class T> struct Counting {
  std::optional<std::vector<T>> tmult;
  std::vector<T> pcoef;
  std::vector<T> pdrop;
  std::string arcsLost;

  static Counting identity(size_t places, size_t transitions) {
    Counting c;
    c.tmult = std::vector<T>(transitions, T(0));
    c.pcoef.assign(places, T(0));
    return c;
  }

  void dropArcs(const char* why) {
    if (tmult) { tmult.reset(); arcsLost = why; }
  }
  /** The survivor stands for the dropped transition too. */
  void fused(size_t dropped, size_t survivor) {
    if (tmult) (*tmult)[survivor] = petri::addExact((*tmult)[survivor], petri::addExact((*tmult)[dropped], T(1)));
  }
  /** A constant place holding `marking` tokens leaves the net. */
  void constantDropped(T marking) {
    if (marking > 0) pdrop.push_back(marking);
  }
  /** `kept` absorbs `other`: their coefficients add, and the moves inside
   * the component are no longer moves of this net. */
  void placesFused(size_t other, size_t kept) {
    pcoef[kept] = petri::addExact(pcoef[kept], petri::addExact(pcoef[other], T(1)));
    dropArcs("free components fused: the moves inside them are not the moves of this net");
  }

  /** The record of the published net: retired slots leave, live ones renumber. */
  Counting compact(const std::vector<size_t>& placeMap, const std::vector<size_t>& transitionMap) const {
    constexpr size_t absent = std::numeric_limits<size_t>::max();
    Counting out;
    out.pdrop = pdrop; out.arcsLost = arcsLost;
    size_t live = 0;
    for (size_t p = 0; p < placeMap.size(); ++p) if (placeMap[p] != absent) ++live;
    out.pcoef.assign(live, T(0));
    for (size_t p = 0; p < placeMap.size(); ++p) if (placeMap[p] != absent) out.pcoef[placeMap[p]] = pcoef[p];
    if (tmult) {
      live = 0;
      for (size_t t = 0; t < transitionMap.size(); ++t) if (transitionMap[t] != absent) ++live;
      out.tmult = std::vector<T>(live, T(0));
      for (size_t t = 0; t < transitionMap.size(); ++t)
        if (transitionMap[t] != absent) (*out.tmult)[transitionMap[t]] = (*tmult)[t];
    }
    return out;
  }
  bool weighted() const { for (T k : pcoef) if (k != 0) return true; return false; }
};
}
