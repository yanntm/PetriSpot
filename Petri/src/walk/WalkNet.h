/*
 * WalkNet.h
 *
 * Read-only compiled view of a SparsePetriNet for explicit exploration:
 * per-transition effects and the place-to-transition consumer index sorted
 * by arc weight. See algorithm.md.
 */
#ifndef PETRI_WALK_WALKNET_H_
#define PETRI_WALK_WALKNET_H_

#include <algorithm>
#include <cstdint>
#include <vector>

#include "core/MatrixCol.h"
#include "core/SparsePetriNet.h"

namespace petri::walk
{

template<typename T>
  class WalkNet
  {
  public:
    struct Consumer
    {
      uint32_t transition;
      T weight;
    };

  private:
    const SparsePetriNet<T> &net;
    MatrixCol<T> effects;
    std::vector<std::vector<Consumer>> consumers;
    std::vector<T> maxConsumeWeight;

  public:
    explicit WalkNet (const SparsePetriNet<T> &n)
        : net (n), effects (n.getPlaceCount (), 0),
          consumers (n.getPlaceCount ()),
          maxConsumeWeight (n.getPlaceCount (), 0)
    {
      const MatrixCol<T> &pt = net.getFlowPT ();
      const MatrixCol<T> &tp = net.getFlowTP ();
      const size_t nbT = net.getTransitionCount ();
      effects.reserveColumns (nbT);
      for (size_t t = 0; t < nbT; ++t) {
        effects.appendColumn (
            SparseArray<T>::sumProd (-1, pt.getColumn (t), 1, tp.getColumn (t)));
        const SparseArray<T> &pre = pt.getColumn (t);
        for (size_t i = 0; i < pre.size (); ++i) {
          size_t p = pre.keyAt (i);
          T w = pre.valueAt (i);
          consumers[p].push_back ({ static_cast<uint32_t> (t), w });
          maxConsumeWeight[p] = std::max (maxConsumeWeight[p], w);
        }
      }
      for (auto &cs : consumers) {
        std::sort (cs.begin (), cs.end (), [] (const Consumer &a, const Consumer &b) {
          return a.weight < b.weight || (a.weight == b.weight && a.transition < b.transition);
        });
      }
    }

    const SparsePetriNet<T>& getNet () const
    {
      return net;
    }
    size_t transitionCount () const
    {
      return net.getTransitionCount ();
    }
    size_t placeCount () const
    {
      return net.getPlaceCount ();
    }
    const SparseArray<T>& pre (size_t t) const
    {
      return net.getFlowPT ().getColumn (t);
    }
    const SparseArray<T>& effect (size_t t) const
    {
      return effects.getColumn (t);
    }
    bool hasEffect (size_t t) const
    {
      return effects.getColumn (t).size () > 0;
    }
    const std::vector<Consumer>& consumersOf (size_t p) const
    {
      return consumers[p];
    }
    T maxWeight (size_t p) const
    {
      return maxConsumeWeight[p];
    }
    const std::vector<T>& initialMarking () const
    {
      return net.getMarks ();
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_WALKNET_H_ */
