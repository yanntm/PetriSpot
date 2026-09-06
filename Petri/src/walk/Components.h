/*
 * Components.h
 *
 * The net as a set of sequential processes, read from its P-flows: a flow
 * with non-negative coefficients is a component, its places the local states,
 * its value at the initial marking the number of tokens it carries (one, for
 * a process). A transition that takes from a component and gives back to it
 * is a local move of that process; a transition touching several components
 * synchronises them. The local graph of a component gives exact local
 * distances (shortest paths) to any of its places, the guide of a quest that
 * drives the process there. See algorithm.md.
 */
#ifndef PETRI_WALK_COMPONENTS_H_
#define PETRI_WALK_COMPONENTS_H_

#include <algorithm>
#include <cstdint>
#include <deque>
#include <limits>
#include <map>
#include <memory>
#include <tuple>
#include <mutex>
#include <ostream>
#include <unordered_map>
#include <vector>

#include "core/MatrixCol.h"
#include "core/SparseArray.h"
#include "walk/WalkNet.h"

namespace petri::walk
{

constexpr uint32_t UNREACHED = std::numeric_limits<uint32_t>::max ();

template<typename T>
  class Components
  {
  public:
    struct Component
    {
      std::vector<size_t> places;       // global place ids, sorted
      std::vector<T> weights;           // the flow's coefficient per place
      T value = 0;                      // weighted token count, invariant along every run
      std::vector<std::vector<std::pair<uint32_t, uint32_t>>> preds; // local q -> (local p, transition): p moves to q by t
      uint64_t signature = 0;           // equal for components whose local graphs look alike (degree profile)
      uint32_t isoClass = 0;

      size_t size () const
      {
        return places.size ();
      }
      /** Local index of a global place, or -1. */
      int32_t indexOf (size_t p) const
      {
        auto it = std::lower_bound (places.begin (), places.end (), p);
        return it != places.end () && *it == p ? static_cast<int32_t> (it - places.begin ()) : -1;
      }
    };

  private:
    std::vector<Component> comps;
    std::vector<std::vector<uint32_t>> ofPlace;   // place -> components containing it, smallest first
    std::vector<uint32_t> syncDegree;             // per transition: components it touches
    mutable std::map<std::tuple<uint32_t, size_t, uint32_t>, std::shared_ptr<const std::vector<uint32_t>>> distCache;
    mutable std::mutex cacheMutex;   // strategies of several threads ask at once

  public:
    /** From the columns of a flow basis: one component per non-negative column. */
    Components (const WalkNet<T> &net, const MatrixCol<T> &flows)
        : ofPlace (net.placeCount ()), syncDegree (net.transitionCount (), 0)
    {
      const std::vector<T> &m0 = net.initialMarking ();
      for (size_t c = 0; c < flows.getColumnCount (); ++c) {
        const SparseArray<T> &col = flows.getColumn (c);
        bool positive = true;
        for (size_t i = 0; i < col.size (); ++i) positive = positive && col.valueAt (i) > 0;
        if (!positive || col.size () == 0) continue;
        Component k;
        for (size_t i = 0; i < col.size (); ++i) {
          k.places.push_back (col.keyAt (i));
          k.weights.push_back (col.valueAt (i));
          k.value += col.valueAt (i) * m0[col.keyAt (i)];
        }
        // keyAt is increasing: places sorted already
        k.preds.assign (k.places.size (), {});
        comps.push_back (std::move (k));
      }
      for (uint32_t c = 0; c < comps.size (); ++c)
        for (size_t p : comps[c].places) ofPlace[p].push_back (c);
      for (auto &l : ofPlace)
        std::sort (l.begin (), l.end (), [this] (uint32_t a, uint32_t b) { return comps[a].size () < comps[b].size (); });
      buildEdges (net);
      sign ();
    }

    size_t size () const
    {
      return comps.size ();
    }
    const Component& component (uint32_t c) const
    {
      return comps[c];
    }
    /** Components containing p, the most specific (smallest) first; empty when p is in no flow. */
    const std::vector<uint32_t>& componentsOf (size_t p) const
    {
      return ofPlace[p];
    }
    uint32_t syncDegreeOf (uint32_t t) const
    {
      return syncDegree[t];
    }

    /**
     * Shortest path, in moves, from every place of component c to target,
     * using only the transitions that synchronise at most maxSync components
     * (1: the moves the process makes alone); UNREACHED where none. Cached: a
     * quest asks for one table per atom.
     */
    std::shared_ptr<const std::vector<uint32_t>> distancesTo (uint32_t c, size_t target, uint32_t maxSync) const
    {
      std::lock_guard<std::mutex> lock (cacheMutex);
      auto key = std::make_tuple (c, target, maxSync);
      auto it = distCache.find (key);
      if (it != distCache.end ()) return it->second;
      const Component &k = comps[c];
      auto dist = std::make_shared<std::vector<uint32_t>> (k.size (), UNREACHED);
      int32_t ti = k.indexOf (target);
      if (ti >= 0) {
        std::deque<uint32_t> queue;
        (*dist)[ti] = 0;
        queue.push_back (static_cast<uint32_t> (ti));
        while (!queue.empty ()) {
          uint32_t q = queue.front ();
          queue.pop_front ();
          for (const auto &e : k.preds[q]) {
            if (syncDegree[e.second] > maxSync) continue;
            if ((*dist)[e.first] == UNREACHED) {
              (*dist)[e.first] = (*dist)[q] + 1;
              queue.push_back (e.first);
            }
          }
        }
      }
      distCache[key] = dist;
      return dist;
    }

    /**
     * The guide of a quest: distances by the process's own moves where they
     * exist, and where a place only reaches the target through
     * synchronisations, the distance over every move plus an offset, so that
     * a process that can walk alone is driven first and one behind a barrier
     * still knows which way to lean.
     */
    std::shared_ptr<const std::vector<uint32_t>> questDistances (uint32_t c, size_t target) const
    {
      auto local = distancesTo (c, target, 1);
      auto full = distancesTo (c, target, std::numeric_limits<uint32_t>::max ());
      auto out = std::make_shared<std::vector<uint32_t>> (*local);
      for (size_t i = 0; i < out->size (); ++i)
        if ((*out)[i] == UNREACHED && (*full)[i] != UNREACHED) (*out)[i] = (*full)[i] + BARRIER_OFFSET;
      return out;
    }

    static constexpr uint32_t BARRIER_OFFSET = 1000;

    /** A few numbers on the decomposition: components, sizes, isomorphism classes, synchronisation degrees. */
    void printStats (std::ostream &os) const
    {
      std::map<size_t, size_t> sizes;
      std::map<uint32_t, size_t> classes;
      for (const auto &k : comps) {
        ++sizes[k.size ()];
        ++classes[k.isoClass];
      }
      std::map<uint32_t, size_t> degrees;
      for (uint32_t d : syncDegree) ++degrees[d];
      os << "Components: " << comps.size () << " (sizes";
      for (const auto &s : sizes) os << " " << s.first << "x" << s.second;
      os << "), " << classes.size () << " isomorphism classes by degree profile; transitions by components synchronised:";
      for (const auto &d : degrees) os << " " << d.first << ":" << d.second;
      os << std::endl;
    }

  private:
    /** A transition taking from a component and giving back to it moves the process: one edge per (from, to) pair. */
    void buildEdges (const WalkNet<T> &net)
    {
      std::unordered_map<uint32_t, std::pair<std::vector<uint32_t>, std::vector<uint32_t>>> touched; // comp -> (from local, to local)
      for (uint32_t t = 0; t < net.transitionCount (); ++t) {
        touched.clear ();
        const SparseArray<T> &eff = net.effect (t);
        for (size_t i = 0; i < eff.size (); ++i) {
          size_t p = eff.keyAt (i);
          for (uint32_t c : ofPlace[p]) {
            int32_t li = comps[c].indexOf (p);
            if (eff.valueAt (i) < 0) touched[c].first.push_back (static_cast<uint32_t> (li));
            else touched[c].second.push_back (static_cast<uint32_t> (li));
          }
        }
        syncDegree[t] = static_cast<uint32_t> (touched.size ());
        for (const auto &ct : touched)
          for (uint32_t from : ct.second.first)
            for (uint32_t to : ct.second.second) comps[ct.first].preds[to].emplace_back (from, t);
      }
    }

    /** Degree profile hash: a cheap filter for isomorphic processes, exact labelling can come later. */
    void sign ()
    {
      std::map<uint64_t, uint32_t> classes;
      for (auto &k : comps) {
        std::vector<uint64_t> profile;
        std::vector<uint32_t> outdeg (k.size (), 0);
        for (uint32_t q = 0; q < k.size (); ++q)
          for (const auto &e : k.preds[q]) ++outdeg[e.first];
        for (uint32_t q = 0; q < k.size (); ++q) profile.push_back ((static_cast<uint64_t> (k.preds[q].size ()) << 32) | outdeg[q]);
        std::sort (profile.begin (), profile.end ());
        uint64_t h = 1469598103934665603ull ^ k.size ();
        for (uint64_t v : profile) h = (h ^ v) * 1099511628211ull;
        k.signature = h;
        auto it = classes.find (h);
        if (it == classes.end ()) it = classes.emplace (h, static_cast<uint32_t> (classes.size ())).first;
        k.isoClass = it->second;
      }
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_COMPONENTS_H_ */
