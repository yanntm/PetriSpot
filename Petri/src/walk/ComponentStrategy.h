/*
 * ComponentStrategy.h
 *
 * A quest per atom, in stages. The goal, a conjunction of `place >= k` atoms,
 * is read as one process per atom to drive to a place. When the places lie
 * behind barriers (transitions synchronising several processes), a stage is
 * set: the barrier every process can reach alone whose crossing brings the
 * goal closest, its pre-places the quests; when it becomes enabled it is
 * fired and the next stage is chosen, and a barrier that did not help is not
 * chosen again in the run. Within a stage the distance
 * of the marking is the sum over the unsatisfied quests of the local shortest
 * path of the quest's process to its place, the choice is greedy on that
 * distance among a few sampled enabled transitions, with an epsilon share of
 * uniform moves, and a satisfied quest freezes its process: a transition that
 * would take the tokens away is refused, so the processes arrive one by one
 * and wait, which turns the product of their arrival chances into a sum.
 * See algorithm.md.
 */
#ifndef PETRI_WALK_COMPONENTSTRATEGY_H_
#define PETRI_WALK_COMPONENTSTRATEGY_H_

#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <ostream>
#include <string>
#include <tuple>
#include <vector>

#include "expr/Expression.h"
#include "walk/Components.h"
#include "walk/Strategy.h"

namespace petri::walk
{

template<typename T>
  class ComponentStrategy : public Strategy<T>
  {
    struct Quest
    {
      size_t place;
      long long need;                     // tokens required on place
      uint32_t comp;                      // the most specific component holding place, NONE when no component does
      const std::vector<uint32_t> *dist;  // quest distances to place by local index, null without a component
    };

    static constexpr uint32_t NONE = std::numeric_limits<uint32_t>::max ();

    const WalkNet<T> &net;
    const Components<T> &comps;
    std::vector<Quest> goalQuests;        // the atoms of the goal
    std::vector<Quest> quests;            // the current stage: the goal's atoms, or a barrier's pre-places
    uint32_t stageBarrier = NONE;         // the barrier the stage prepares, NONE when the stage is the goal
    bool needReplan = true;
    uint64_t hBefore = 0;                 // goal distance when the current barrier was chosen
    std::vector<uint32_t> tabu;           // barriers whose crossing did not lower the goal distance, this run
    unsigned epsilonPercent;
    size_t sampleSize;
    uint64_t stallLimit;
    uint64_t sinceImprovement = 0;

    // per step scratch
    std::vector<uint64_t> cur;            // current distance per quest, 0 when satisfied
    std::vector<bool> frozenPlace;        // per place: a satisfied quest holds it
    std::vector<size_t> frozenList;

  public:
    /** Print the quests and their distances from the initial marking on stderr. */
    void describe (std::ostream &os, const Marking<T> &m0) const
    {
      os << "sync: " << goalQuests.size () << " quests" << std::endl;
      for (const Quest &q : goalQuests) {
        if (q.comp == NONE) {
          os << "  place " << q.place << " >= " << q.need << " in no component" << std::endl;
          continue;
        }
        const auto &k = comps.component (q.comp);
        uint64_t d = UNREACHED;
        for (size_t j = 0; j < k.size (); ++j)
          if (m0.get (k.places[j]) > 0 && (*q.dist)[j] < d) d = (*q.dist)[j];
        os << "  place " << q.place << " >= " << q.need << " component " << q.comp << " (" << k.size ()
            << " places, " << k.value << " tokens) distance from the initial marking "
            << (d == UNREACHED ? std::string ("unreached") : std::to_string (d)) << std::endl;
      }
    }

    /** Print the stage on stderr at the onset of the first debugSteps stalls. */
    uint64_t debugSteps = 0;

    /** The stage as it stands: the barrier, each quest with its distance and the marked places of its process. */
    void describeState (std::ostream &os, const Marking<T> &m, const EnabledSet<T> &enabled) const
    {
      os << "sync stage: barrier " << (stageBarrier == NONE ? std::string ("goal") : std::to_string (stageBarrier))
          << ", " << enabled.size () << " enabled, refusals " << refusals << std::endl;
      for (size_t i = 0; i < quests.size (); ++i) {
        const Quest &q = quests[i];
        if (q.comp == NONE) {
          os << "  place " << q.place << " >= " << q.need << " in no component, cur " << cur[i] << std::endl;
          continue;
        }
        const auto &k = comps.component (q.comp);
        os << "  place " << q.place << " >= " << q.need << " component " << q.comp << " cur " << cur[i] << " marked:";
        for (size_t j = 0; j < k.size (); ++j)
          if (m.get (k.places[j]) > 0) os << " p" << k.places[j] << "(d" << (*q.dist)[j] << ")";
        os << std::endl;
        if (cur[i] == 0) continue;
        // why the process does not move: the consumers of its marked places
        for (size_t j = 0; j < k.size (); ++j) {
          if (m.get (k.places[j]) == 0) continue;
          const auto &cs = net.consumersOf (k.places[j]);
          os << "    consumers of p" << k.places[j] << " (" << cs.size () << "):";
          size_t shown = 0;
          for (const auto &c : cs) {
            if (shown++ >= 6) { os << " ..."; break; }
            uint32_t t = c.transition;
            const SparseArray<T> &eff = net.effect (t);
            uint64_t after = UNREACHED;
            for (size_t e = 0; e < eff.size (); ++e) {
              if (eff.valueAt (e) <= 0) continue;
              int32_t li = k.indexOf (eff.keyAt (e));
              if (li >= 0 && (*q.dist)[li] < after) after = (*q.dist)[li];
            }
            os << " t" << t << "[pre " << net.pre (t).size () << ", sync " << comps.syncDegreeOf (t)
                << (enabled.isEnabled (t) ? ", enabled" : ", unsat " + std::to_string (enabled.unsatOf (t)))
                << (breaksFrozen (t) ? ", refused" : "") << ", to d" << (after == UNREACHED ? std::string ("-") : std::to_string (after)) << "]";
          }
          os << std::endl;
        }
      }
    }

    /** Instrumentation. */
    uint64_t minDistance = std::numeric_limits<uint64_t>::max ();
    uint64_t minDistanceThisRun = std::numeric_limits<uint64_t>::max ();
    uint64_t refusals = 0;                // candidates refused for breaking a frozen process
    uint64_t replans = 0;                 // stages set up
    uint64_t barriersFired = 0;           // stage barriers fired
    uint64_t hopeless = 0;                // runs abandoned with a process where its place is unreachable
    SparseArray<T> bestMarkingThisRun;

    /**
     * Reads the goal; supported() says whether every atom was a `place >= k`
     * with a component. The caller falls back to another strategy otherwise.
     */
    ComponentStrategy (const WalkNet<T> &n, const Components<T> &c, const petri::expr::Expression &goal,
                       unsigned epsilon, size_t sample, uint64_t stall)
        : net (n), comps (c), epsilonPercent (epsilon), sampleSize (sample == 0 ? 16 : sample),
          stallLimit (stall), frozenPlace (n.placeCount (), false)
    {
      ok = collect (goal);
      goalQuests = quests;
      cur.assign (quests.size (), 0);
    }

    bool supported () const
    {
      return ok && !quests.empty ();
    }
    size_t questCount () const
    {
      return quests.size ();
    }

    void onReset () override
    {
      minDistanceThisRun = std::numeric_limits<uint64_t>::max ();
      sinceImprovement = 0;
      needReplan = true;
      tabu.clear ();
    }

    bool bestOfRun (SparseArray<T> &m, uint64_t &h) const override
    {
      if (minDistanceThisRun == std::numeric_limits<uint64_t>::max ()) return false;
      m = bestMarkingThisRun;
      h = minDistanceThisRun;
      return true;
    }

    uint32_t choose (WalkContext<T> &ctx) override
    {
      size_t n = ctx.enabled.size ();
      if (stallLimit && sinceImprovement >= stallLimit) return RESTART;
      if (needReplan && !replan (ctx.marking)) return RESTART; // no barrier and no goal within the processes' reach
      if (stageBarrier != NONE && ctx.enabled.isEnabled (stageBarrier)) {
        // the processes are aligned: cross the barrier, then plan the next stage
        ++barriersFired;
        needReplan = true;
        return stageBarrier;
      }
      uint64_t total = measure (ctx.marking);
      if (total >= HOPELESS) {
        // a process sits where it can never reach its place: this run is over
        ++hopeless;
        return RESTART;
      }
      if (total < minDistanceThisRun) {
        minDistanceThisRun = total;
        sinceImprovement = 0;
        bestMarkingThisRun = ctx.marking.sparse ();
        if (total < minDistance) minDistance = total;
      } else {
        ++sinceImprovement;
      }
      bool trace = debugSteps > 0 && total <= 2;
      if (trace) {
        --debugSteps;
        describeState (std::cerr, ctx.marking, ctx.enabled);
      }
      size_t k = sampleSize >= n ? n : sampleSize;
      if (n == 1 || (epsilonPercent > 0 && ctx.rng () % 100 < epsilonPercent)) {
        // a random move, but never one that undoes a process already in place
        for (size_t i = 0; i < k; ++i) {
          uint32_t t = ctx.enabled.at (ctx.rng () % n);
          if (!breaksFrozen (t)) return t;
        }
        return ctx.enabled.at (ctx.rng () % n);
      }

      uint32_t best = RESTART;
      uint64_t bestScore = std::numeric_limits<uint64_t>::max ();
      uint32_t ties = 0;
      for (size_t i = 0; i < k; ++i) {
        uint32_t t = k == n ? ctx.enabled.at (i) : ctx.enabled.at (ctx.rng () % n);
        if (breaksFrozen (t)) {
          ++refusals;
          continue;
        }
        uint64_t s = score (t, total);
        if (s < bestScore) {
          bestScore = s;
          best = t;
          ties = 1;
        } else if (s == bestScore && ctx.rng () % ++ties == 0) {
          best = t;
        }
      }
      if (trace) {
        std::cerr << "  chose " << (best == RESTART ? std::string ("fallback") : "t" + std::to_string (best))
            << " score " << bestScore << " of " << k << " sampled; enabled:";
        for (size_t i = 0; i < std::min<size_t> (n, 12); ++i) {
          uint32_t t = ctx.enabled.at (i);
          std::cerr << " t" << t << (breaksFrozen (t) ? "(refused)" : "") << "=" << score (t, total);
        }
        std::cerr << std::endl;
      }
      // every candidate would undo a process that has arrived: move anyway rather than stand still
      if (best == RESTART) return ctx.enabled.at (ctx.rng () % n);
      return best;
    }

  private:
    bool ok = true;
    static constexpr uint64_t HOPELESS = 1000000; // distance of a quest whose process cannot reach its place

    /**
     * Set the stage from m. The goal's own atoms when every process can reach
     * its place alone. Otherwise the barrier (a transition synchronising
     * several processes) that every process can reach alone and whose firing
     * brings the goal closest: the summed component distance of the marking it
     * produces to the goal places, with the offset for a process still behind
     * another barrier; the cheapest to align among equals. A barrier whose
     * crossing did not lower that sum is not chosen again this run. False when
     * no barrier can be reached and the goal cannot either.
     */
    bool replan (const Marking<T> &m)
    {
      needReplan = false;
      ++replans;
      locate (m);
      // where the goal stands from here, per quest
      uint64_t goalNow = 0;
      bool goalLocal = true;
      goalCur.resize (goalQuests.size ());
      for (size_t i = 0; i < goalQuests.size (); ++i) {
        goalCur[i] = distanceOf (m, goalQuests[i]);
        if (goalCur[i] >= Components<T>::BARRIER_OFFSET) goalLocal = false;
        goalNow += goalCur[i];
      }
      if (stageBarrier != NONE && goalNow >= hBefore) tabu.push_back (stageBarrier);
      hBefore = goalNow;
      uint32_t best = NONE;
      uint64_t bestAfter = std::numeric_limits<uint64_t>::max ();
      uint64_t bestAlign = std::numeric_limits<uint64_t>::max ();
      if (!goalLocal) {
        for (uint32_t t = 0; t < net.transitionCount (); ++t) {
          if (comps.syncDegreeOf (t) <= 1) continue;
          if (std::find (tabu.begin (), tabu.end (), t) != tabu.end ()) continue;
          uint64_t align = alignCost (t);
          if (align == UNREACHED) continue;
          uint64_t after = goalAfter (t);
          if (std::tie (after, align) < std::tie (bestAfter, bestAlign)) {
            bestAfter = after;
            bestAlign = align;
            best = t;
          }
        }
        if (best == NONE) return false;
      }
      stageBarrier = best;
      if (debugSteps > 0)
        std::cerr << "replan goal " << goalNow << (goalLocal ? " local" : "") << ", stage "
            << (best == NONE ? std::string ("goal") : "t" + std::to_string (best)) << " after " << bestAfter
            << " align " << bestAlign << " tabu " << tabu.size () << std::endl;
      quests.clear ();
      if (best == NONE) {
        quests = goalQuests;
      } else {
        const SparseArray<T> &pre = net.pre (best);
        for (size_t i = 0; i < pre.size (); ++i) addQuest (pre.keyAt (i), pre.valueAt (i));
      }
      cur.assign (quests.size (), 0);
      return true;
    }

    /** Distance of quest q from m: 0 when satisfied, the process's distance to the place plus one otherwise. */
    uint64_t distanceOf (const Marking<T> &m, const Quest &q) const
    {
      if (m.get (q.place) >= q.need) return 0;
      if (q.comp == NONE) return 1; // unguided: a move away, as far as we know
      const auto &k = comps.component (q.comp);
      uint64_t d = UNREACHED;
      for (size_t j = 0; j < k.size (); ++j)
        if (m.get (k.places[j]) > 0 && (*q.dist)[j] < d) d = (*q.dist)[j];
      return d == UNREACHED ? HOPELESS : d + 1;
    }

    /**
     * Summed goal distance of the marking t produces from the located marking:
     * a process t moves stands at t's post place in its component, the others
     * where they are (goalNow holds their distances, computed by replan).
     */
    uint64_t goalAfter (uint32_t t) const
    {
      uint64_t s = 0;
      for (size_t i = 0; i < goalQuests.size (); ++i) {
        const Quest &q = goalQuests[i];
        uint64_t after = UNREACHED;
        bool moved = false;
        if (q.comp != NONE) {
          for (const auto &cl : comps.producedBy (t)) {
            if (cl.first != q.comp) continue;
            moved = true;
            if ((*q.dist)[cl.second] < after) after = (*q.dist)[cl.second];
          }
        }
        if (!moved) s += goalCur[i];
        else s += after == UNREACHED ? HOPELESS : (after == 0 ? 0 : after + 1);
      }
      return s;
    }

    // Where each process stands, refreshed once per stage choice: the marked
    // local places of every component, so that a distance is one table lookup.
    std::vector<std::vector<uint32_t>> standing;
    std::vector<uint64_t> goalCur;         // the goal quests' distances at the stage choice

    void locate (const Marking<T> &m)
    {
      standing.assign (comps.size (), {});
      for (uint32_t c = 0; c < comps.size (); ++c) {
        const auto &k = comps.component (c);
        for (uint32_t j = 0; j < k.size (); ++j)
          if (m.get (k.places[j]) > 0) standing[c].push_back (j);
      }
    }

    /** Distance of the process holding p to p from where it stands (after locate); UNREACHED when no component holds p. */
    uint64_t localDistance (size_t p) const
    {
      const auto &cs = comps.componentsOf (p);
      if (cs.empty ()) return UNREACHED;
      const std::vector<uint32_t> &dist = comps.questDistances (cs[0], p);
      uint64_t d = UNREACHED;
      for (uint32_t j : standing[cs[0]])
        if (dist[j] < d) d = dist[j];
      return d;
    }

    /** The cost of aligning every process on the pre-places of t, UNREACHED when one cannot get there alone. */
    uint64_t alignCost (uint32_t t) const
    {
      const SparseArray<T> &pre = net.pre (t);
      uint64_t align = 0;
      for (size_t i = 0; i < pre.size (); ++i) {
        uint64_t d = localDistance (pre.keyAt (i));
        if (d >= Components<T>::BARRIER_OFFSET) return UNREACHED;
        align += d;
      }
      return align;
    }

    /** A quest on place >= need; without a component to guide it the quest only says satisfied or not. */
    void addQuest (size_t p, long long need)
    {
      const auto &cs = comps.componentsOf (p);
      if (cs.empty ()) quests.push_back ({ p, need, NONE, nullptr });
      else quests.push_back ({ p, need, cs[0], &comps.questDistances (cs[0], p) });
    }

    /** Every atom of a conjunction as a quest; anything else is unsupported. */
    bool collect (const petri::expr::Expression &e)
    {
      using petri::expr::Expression;
      using petri::expr::Cmp;
      switch (e.kind) {
      case Expression::Kind::And:
        for (const auto &c : e.children)
          if (!collect (c)) return false;
        return true;
      case Expression::Kind::Atom: {
        const auto &a = e.atom;
        if (a.terms.size () != 1 || a.terms[0].second <= 0) return false;
        long long need;
        if (a.op == Cmp::GE) need = (a.constant + a.terms[0].second - 1) / a.terms[0].second;
        else if (a.op == Cmp::GT) need = a.constant / a.terms[0].second + 1;
        else if (a.op == Cmp::EQ) need = a.constant / a.terms[0].second;
        else return false;
        addQuest (a.terms[0].first, need);
        return true;
      }
      default:
        return false;
      }
    }

    /** Sum of the local distances of the unsatisfied quests; refreshes cur and the frozen places. */
    uint64_t measure (const Marking<T> &m)
    {
      for (size_t p : frozenList) frozenPlace[p] = false;
      frozenList.clear ();
      uint64_t total = 0;
      for (size_t i = 0; i < quests.size (); ++i) {
        const Quest &q = quests[i];
        cur[i] = distanceOf (m, q);
        if (cur[i] == 0) {
          frozenPlace[q.place] = true;
          frozenList.push_back (q.place);
        }
        total += cur[i];
      }
      return total;
    }

    /** Would t take tokens from a place a satisfied quest holds. */
    bool breaksFrozen (uint32_t t) const
    {
      if (frozenList.empty ()) return false;
      const SparseArray<T> &eff = net.effect (t);
      for (size_t i = 0; i < eff.size (); ++i)
        if (eff.valueAt (i) < 0 && frozenPlace[eff.keyAt (i)]) return true;
      return false;
    }

    /** Distance after firing t: only the quests whose component t feeds change. */
    uint64_t score (uint32_t t, uint64_t total) const
    {
      const SparseArray<T> &eff = net.effect (t);
      uint64_t s = total;
      for (size_t i = 0; i < quests.size (); ++i) {
        if (cur[i] == 0) continue;
        const Quest &q = quests[i];
        if (q.comp == NONE) {
          // unguided: a transition feeding the place is the only clue
          if (net.effect (t).get (q.place) > 0) s -= 1;
          continue;
        }
        const auto &k = comps.component (q.comp);
        uint64_t after = UNREACHED;
        bool feeds = false;
        for (size_t j = 0; j < eff.size (); ++j) {
          if (eff.valueAt (j) <= 0) continue;
          int32_t li = k.indexOf (eff.keyAt (j));
          if (li < 0) continue;
          feeds = true;
          if ((*q.dist)[li] < after) after = (*q.dist)[li];
        }
        if (!feeds) continue;               // t does not move this process
        // a move to a place from which the target is unreachable strands the process
        uint64_t nd = after == UNREACHED ? HOPELESS : after + 1;
        if (k.value > 1 && cur[i] < nd) nd = cur[i]; // several tokens: the nearest still counts
        s = s - cur[i] + nd;
      }
      return s;
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_COMPONENTSTRATEGY_H_ */
