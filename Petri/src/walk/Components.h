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
      std::vector<std::vector<uint32_t>> quest; // [local target][local from]: the quest distance (see questDistances)
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
    // per transition: the places it consumes from and produces into, as (component, local index),
    // one entry per component the place belongs to; the tables a stage choice scans
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> consumed, produced;

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
        // a single place is a constant, not a process
        if (!positive || col.size () < 2) continue;
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
      for (uint32_t c = 0; c < comps.size (); ++c) {
        Component &k = comps[c];
        k.quest.resize (k.size ());
        for (uint32_t target = 0; target < k.size (); ++target) k.quest[target] = questDistancesOf (k, target);
      }
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
    /** The (component, local place) pairs t takes tokens from. */
    const std::vector<std::pair<uint32_t, uint32_t>>& consumedBy (uint32_t t) const
    {
      return consumed[t];
    }
    /** The (component, local place) pairs t puts tokens into. */
    const std::vector<std::pair<uint32_t, uint32_t>>& producedBy (uint32_t t) const
    {
      return produced[t];
    }

    /**
     * The quest distances of component c to its local place target, indexed by
     * the local place a process stands on: its shortest path by its own moves
     * (transitions synchronising no other component) where one exists, and
     * where only synchronisations lead there, the distance over every move
     * plus BARRIER_OFFSET, so that a process that can walk alone is driven
     * first and one behind a barrier still knows which way to lean. UNREACHED
     * where no path at all exists. Computed once for every (component, place).
     */
    const std::vector<uint32_t>& questDistances (uint32_t c, size_t target) const
    {
      return comps[c].quest[static_cast<size_t> (comps[c].indexOf (target))];
    }

    static constexpr uint32_t BARRIER_OFFSET = 1000;

  private:
    /** Breadth first over the reversed local edges, keeping the transitions of sync degree at most maxSync. */
    static std::vector<uint32_t> distancesTo (const Component &k, uint32_t target, uint32_t maxSync,
                                              const std::vector<uint32_t> &syncDegree)
    {
      std::vector<uint32_t> dist (k.size (), UNREACHED);
      std::deque<uint32_t> queue;
      dist[target] = 0;
      queue.push_back (target);
      while (!queue.empty ()) {
        uint32_t q = queue.front ();
        queue.pop_front ();
        for (const auto &e : k.preds[q]) {
          if (syncDegree[e.second] > maxSync) continue;
          if (dist[e.first] == UNREACHED) {
            dist[e.first] = dist[q] + 1;
            queue.push_back (e.first);
          }
        }
      }
      return dist;
    }

    std::vector<uint32_t> questDistancesOf (const Component &k, uint32_t target) const
    {
      std::vector<uint32_t> local = distancesTo (k, target, 1, syncDegree);
      std::vector<uint32_t> full = distancesTo (k, target, std::numeric_limits<uint32_t>::max (), syncDegree);
      for (size_t i = 0; i < local.size (); ++i)
        if (local[i] == UNREACHED && full[i] != UNREACHED) local[i] = full[i] + BARRIER_OFFSET;
      return local;
    }

  public:

    /** Transitions synchronising at least two components: the barriers a quest must stage through. */
    size_t barrierCount () const
    {
      size_t n = 0;
      for (uint32_t d : syncDegree) n += d > 1;
      return n;
    }

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
      size_t uncovered = 0;
      for (const auto &l : ofPlace) uncovered += l.empty ();
      os << "Components: " << comps.size () << " (sizes";
      for (const auto &s : sizes) os << " " << s.first << "x" << s.second;
      os << "), " << classes.size () << " isomorphism classes by degree profile, " << uncovered
          << " places in no component; transitions by components synchronised:";
      for (const auto &d : degrees) os << " " << d.first << ":" << d.second;
      os << std::endl;
    }

  private:
    /** A transition taking from a component and giving back to it moves the process: one edge per (from, to) pair. */
    void buildEdges (const WalkNet<T> &net)
    {
      consumed.resize (net.transitionCount ());
      produced.resize (net.transitionCount ());
      std::unordered_map<uint32_t, std::pair<std::vector<uint32_t>, std::vector<uint32_t>>> touched; // comp -> (from local, to local)
      for (uint32_t t = 0; t < net.transitionCount (); ++t) {
        touched.clear ();
        const SparseArray<T> &eff = net.effect (t);
        for (size_t i = 0; i < eff.size (); ++i) {
          size_t p = eff.keyAt (i);
          for (uint32_t c : ofPlace[p]) {
            uint32_t li = static_cast<uint32_t> (comps[c].indexOf (p));
            if (eff.valueAt (i) < 0) {
              touched[c].first.push_back (li);
              consumed[t].emplace_back (c, li);
            } else {
              touched[c].second.push_back (li);
              produced[t].emplace_back (c, li);
            }
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
