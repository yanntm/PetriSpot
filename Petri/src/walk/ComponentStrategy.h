/*
 * ComponentStrategy.h
 *
 * Quests in stages, recursively. The goal, a conjunction of `place >= k`
 * atoms, is one quest per atom: the processes holding the place are driven to
 * put k tokens on it, and a token already there is frozen. Within a stage the
 * distance of the marking is the sum of the quests' distances, the choice is
 * greedy on the distance after firing among a few sampled enabled
 * transitions, with an epsilon share of sideways moves. When a stage stalls,
 * some quest behind a barrier and no progress for a while, a sub-stage is
 * pushed: the barrier every process can reach alone whose outcome brings the
 * stalled stage's quests closest, its pre-places the new quests; it is fired
 * the moment it is enabled and its stage popped; a stage that cannot be
 * staged further pops with its barrier tabu at the parent level; the root
 * failing ends the run. See algorithm.md.
 */
#ifndef PETRI_WALK_COMPONENTSTRATEGY_H_
#define PETRI_WALK_COMPONENTSTRATEGY_H_

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
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
    static constexpr uint32_t NONE = std::numeric_limits<uint32_t>::max ();
    static constexpr uint64_t HOPELESS = 1000000;  // distance of a quest whose tokens cannot reach the place
    static constexpr uint64_t STALL_STEPS = 50;    // steps without a quest crossing its barrier before a stage is staged further
    static constexpr size_t MAX_DEPTH = 8;         // stages stacked at most

    struct Quest
    {
      size_t place;
      long long need;                     // tokens required on place
      uint32_t comp;                      // the most specific component holding place, NONE when no component does
    };

    /** One level: a barrier to enable (NONE at the root, the goal) through its quests. */
    struct Stage
    {
      uint32_t barrier;
      std::vector<Quest> quests;
      std::vector<uint64_t> cur;          // distance per quest at the last measure
      uint64_t best = std::numeric_limits<uint64_t>::max ();
      size_t bestBehind = std::numeric_limits<size_t>::max (); // fewest quests behind a barrier seen
      uint64_t sinceImprovement = 0;      // steps since bestBehind fell (or best, when none is behind)
      uint64_t parentBefore = 0;          // the parent's distance when this stage was pushed
      std::vector<uint32_t> tabu;         // sub-barriers that did not help this stage, this run
    };

    const WalkNet<T> &net;
    const Components<T> &comps;
    std::vector<Quest> goalQuests;
    std::vector<Stage> stack;
    bool pendingPop = false;              // the top barrier was returned to fire: pop at the next choice
    uint32_t justFired = NONE;            // the barrier popped, to judge whether it helped its parent
    uint64_t parentBefore = 0;
    unsigned epsilonPercent;
    size_t sampleSize;
    uint64_t stallLimit;
    uint64_t sinceImprovement = 0;        // of the root, for the global stall

    // per step scratch
    std::vector<long long> frozenNeed;    // per place: tokens a quest wants kept there, 0 when none
    std::vector<size_t> frozenList;
    std::vector<std::vector<uint32_t>> standing; // per component: the marked local places, at a stage choice
    mutable std::vector<uint64_t> scratch;
    const Marking<T> *current = nullptr;
    bool ok = true;

  public:
    /**
     * The tokens a single-place atom asks for: `c*p >= k`, `c*p > k` or `c*p == k`
     * read as p >= need (an equality is quested as its lower side; the claim
     * itself is checked on the exact predicate). False for any other atom.
     */
    static bool questNeed (const petri::expr::LinearAtom &a, long long &need)
    {
      using petri::expr::Cmp;
      if (a.terms.size () != 1 || a.terms[0].second <= 0) return false;
      long long c = a.terms[0].second;
      if (a.op == Cmp::GE) need = (a.constant + c - 1) / c;
      else if (a.op == Cmp::GT) need = a.constant / c + 1;
      else if (a.op == Cmp::EQ) need = a.constant / c;
      else return false;
      return true;
    }

    /** Print the stage sequence and the stuck states on stderr, this many times. */
    uint64_t debugSteps = 0;

    /** Instrumentation. */
    uint64_t minDistance = std::numeric_limits<uint64_t>::max ();
    uint64_t minDistanceThisRun = std::numeric_limits<uint64_t>::max ();
    uint64_t refusals = 0;                // candidates refused for breaking a frozen place
    uint64_t replans = 0;                 // stages pushed
    uint64_t barriersFired = 0;           // stage barriers fired
    uint64_t hopeless = 0;                // runs abandoned with a quest that cannot be met
    uint64_t unstageable = 0;             // runs abandoned with the root behind barriers and no barrier to stage
    uint64_t fallbacks = 0;               // stages popped for want of a sub-barrier
    uint64_t maxDepth = 0;                // deepest stack seen
    SparseArray<T> bestMarkingThisRun;

    ComponentStrategy (const WalkNet<T> &n, const Components<T> &c, const petri::expr::Expression &goal,
                       unsigned epsilon, size_t sample, uint64_t stall)
        : net (n), comps (c), epsilonPercent (epsilon), sampleSize (sample == 0 ? 16 : sample), stallLimit (stall),
          frozenNeed (n.placeCount (), 0)
    {
      ok = collect (goal, goalQuests);
      stack.push_back (Stage { NONE, goalQuests, {} });
    }

    bool supported () const
    {
      return ok && !goalQuests.empty ();
    }
    size_t questCount () const
    {
      return goalQuests.size ();
    }

    void onReset () override
    {
      minDistanceThisRun = std::numeric_limits<uint64_t>::max ();
      sinceImprovement = 0;
      stack.resize (1);
      stack[0].best = std::numeric_limits<uint64_t>::max ();
      stack[0].bestBehind = std::numeric_limits<size_t>::max ();
      stack[0].sinceImprovement = 0;
      stack[0].tabu.clear ();
      pendingPop = false;
      justFired = NONE;
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
      current = &ctx.marking;
      if (pendingPop) {
        // the barrier fired (or could not): back to the stage it served
        pendingPop = false;
        if (stack.size () > 1) {
          justFired = stack.back ().barrier;
          parentBefore = stack.back ().parentBefore;
          stack.pop_back ();
        }
      }
      // the top barrier is ready: cross it
      if (stack.back ().barrier != NONE && ctx.enabled.isEnabled (stack.back ().barrier)) {
        ++barriersFired;
        pendingPop = true;
        return stack.back ().barrier;
      }
      // settle the stage to work on: measure, and stage further when stalled
      uint64_t total;
      for (;;) {
        Stage &top = stack.back ();
        total = measure (ctx.marking, top);
        if (total >= HOPELESS) {
          ++hopeless;
          return RESTART;
        }
        size_t open = 0, behind = 0;
        for (uint64_t d : top.cur) {
          open += d > 0;
          behind += d >= Components<T>::BARRIER_OFFSET;
        }
        if (justFired != NONE) {
          // the sub-barrier crossed: did it bring this stage closer? never choose it again here otherwise
          if (total >= parentBefore) top.tabu.push_back (justFired);
          justFired = NONE;
          top.sinceImprovement = 0;
        }
        // progress is a quest crossing its barrier; below that, the distance itself. The far side
        // fluctuates with every local move and would never let a stall be seen.
        if (behind < top.bestBehind || (behind == 0 && total < top.best)) {
          top.bestBehind = behind;
          top.best = std::min (top.best, total);
          top.sinceImprovement = 0;
        } else {
          ++top.sinceImprovement;
        }
        if (stack.size () == 1) {
          if (total < minDistanceThisRun) {
            minDistanceThisRun = total;
            sinceImprovement = 0;
            bestMarkingThisRun = ctx.marking.sparse ();
            if (total < minDistance) minDistance = total;
          } else {
            ++sinceImprovement;
          }
        }
        // a quest behind a barrier cannot be helped by this stage's own moves: when every open
        // quest is, staging is the only move left; when some can still walk, wait for a stall
        if (behind == 0 || (behind < open && top.sinceImprovement < STALL_STEPS)) break;
        // a stage for this stage, or give this stage up
        if (stack.size () < MAX_DEPTH && push (ctx.marking)) continue;
        if (stack.size () == 1) {
          ++unstageable;
          return RESTART;
        }
        ++fallbacks;
        uint32_t failed = top.barrier;
        stack.pop_back ();
        stack.back ().tabu.push_back (failed);
        stack.back ().sinceImprovement = 0;
      }
      Stage &top = stack.back ();
      // a stage pushed just now may already be ready: its barrier is what we want, whatever it consumes
      if (top.barrier != NONE && ctx.enabled.isEnabled (top.barrier)) {
        ++barriersFired;
        pendingPop = true;
        return top.barrier;
      }

      size_t k = sampleSize >= n ? n : sampleSize;
      bool trace = debugSteps > 0 && (total <= 2 || (n <= 2 && stack.size () == 1));
      if (trace) {
        --debugSteps;
        describeState (std::cerr, ctx.marking, ctx.enabled);
      }
      if (n == 1 || (epsilonPercent > 0 && ctx.rng () % 100 < epsilonPercent)) {
        // a sideways move: random among those that neither undo a frozen place nor push a process away
        for (size_t i = 0; i < k; ++i) {
          uint32_t t = ctx.enabled.at (ctx.rng () % n);
          if (!breaksFrozen (t) && score (t, top, total) <= total + 1) return t;
        }
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
        uint64_t s = score (t, top, total);
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
            << " score " << bestScore << " of " << k << " sampled" << std::endl;
      }
      // every candidate would undo a frozen place: move anyway rather than stand still
      if (best == RESTART) return ctx.enabled.at (ctx.rng () % n);
      return best;
    }

    /** The goal's quests and their distances from a marking, on stderr. */
    void describe (std::ostream &os, const Marking<T> &m0) const
    {
      os << "sync: " << goalQuests.size () << " quests" << std::endl;
      for (const Quest &q : goalQuests) {
        if (q.comp == NONE) {
          os << "  place " << q.place << " >= " << q.need << " in no component" << std::endl;
          continue;
        }
        const auto &k = comps.component (q.comp);
        uint64_t d = distanceOf (m0, q);
        os << "  place " << q.place << " >= " << q.need << " component " << q.comp << " (" << k.size ()
            << " places, " << k.value << " tokens) distance from the initial marking "
            << (d >= HOPELESS ? std::string ("unreached") : std::to_string (d)) << std::endl;
      }
    }

    /** The stack of stages as it stands, the top's quests with the marked places of their processes. */
    void describeState (std::ostream &os, const Marking<T> &m, const EnabledSet<T> &enabled) const
    {
      os << "sync stack:";
      for (const Stage &s : stack) os << " " << (s.barrier == NONE ? std::string ("goal") : "t" + std::to_string (s.barrier));
      os << ", " << enabled.size () << " enabled, refusals " << refusals << ", frozen:";
      for (size_t p : frozenList) os << " p" << p << "(keep " << frozenNeed[p] << " of " << m.get (p) << ")";
      os << std::endl;
      const Stage &tp = stack.back ();
      if (tp.barrier != NONE) {
        const SparseArray<T> &pre = net.pre (tp.barrier);
        os << "  barrier t" << tp.barrier << (enabled.isEnabled (tp.barrier) ? " enabled" : " unsat " + std::to_string (enabled.unsatOf (tp.barrier))) << ", pre:";
        for (size_t i = 0; i < pre.size (); ++i) os << " p" << pre.keyAt (i) << ">=" << pre.valueAt (i) << "(" << m.get (pre.keyAt (i)) << ")";
        os << std::endl;
      }
      const Stage &top = stack.back ();
      for (size_t i = 0; i < top.quests.size (); ++i) {
        const Quest &q = top.quests[i];
        os << "  place " << q.place << " >= " << q.need << " cur " << top.cur[i];
        if (q.comp != NONE) {
          os << " marked:";
          for (uint32_t c : comps.componentsOf (q.place)) {
            const auto &k = comps.component (c);
            const std::vector<uint32_t> &dist = comps.questDistances (c, q.place);
            for (size_t j = 0; j < k.size (); ++j)
              if (m.get (k.places[j]) > 0) os << " p" << k.places[j] << "(d" << dist[j] << ")";
          }
        }
        os << std::endl;
      }
    }

  private:
    /** Every atom of a conjunction as a quest; anything else is unsupported. */
    bool collect (const petri::expr::Expression &e, std::vector<Quest> &out)
    {
      using petri::expr::Expression;
      using petri::expr::Cmp;
      switch (e.kind) {
      case Expression::Kind::And:
        for (const auto &c : e.children)
          if (!collect (c, out)) return false;
        return true;
      case Expression::Kind::Atom: {
        long long need;
        if (!questNeed (e.atom, need)) return false;
        addQuest (e.atom.terms[0].first, need, out);
        return true;
      }
      default:
        return false;
      }
    }

    /** A quest on place >= need; without a component to guide it the quest only says satisfied or not. */
    void addQuest (size_t p, long long need, std::vector<Quest> &out) const
    {
      const auto &cs = comps.componentsOf (p);
      out.push_back ({ p, need, cs.empty () ? NONE : cs[0] });
    }

    /**
     * Distance of quest q from m (plus a virtual effect): 0 when satisfied;
     * otherwise the tokens still missing come from the components holding the
     * place, and the distance is the sum of the local paths of the nearest
     * ones, a token per unit of marking, plus one each. Gathering k tokens is
     * a synchronisation of those components with the place.
     */
    uint64_t distanceOf (const Marking<T> &m, const Quest &q, const SparseArray<T> *delta = nullptr) const
    {
      auto at = [&] (size_t place) -> long long {
        long long v = m.get (place);
        if (delta) v += delta->get (place);
        return v;
      };
      long long have = at (q.place);
      if (have >= q.need) return 0;
      if (q.comp == NONE) return 1; // unguided: a move away, as far as we know
      size_t missing = static_cast<size_t> (q.need - have);
      scratch.clear ();
      for (uint32_t c : comps.componentsOf (q.place)) {
        const auto &k = comps.component (c);
        const std::vector<uint32_t> &dist = comps.questDistances (c, q.place);
        for (size_t j = 0; j < k.size (); ++j) {
          if (k.places[j] == q.place || dist[j] == UNREACHED) continue;
          long long v = at (k.places[j]);
          for (long long u = 0; u < v; ++u) scratch.push_back (dist[j]);
        }
      }
      if (scratch.size () < missing) return HOPELESS;
      std::nth_element (scratch.begin (), scratch.begin () + static_cast<long> (missing) - 1, scratch.end ());
      uint64_t sum = 0;
      for (size_t i = 0; i < missing; ++i) sum += scratch[i] + 1;
      return sum;
    }

    /** Sum of the stage's quest distances; refreshes cur and the frozen places. */
    uint64_t measure (const Marking<T> &m, Stage &stage)
    {
      for (size_t p : frozenList) frozenNeed[p] = 0;
      frozenList.clear ();
      stage.cur.resize (stage.quests.size ());
      uint64_t total = 0;
      for (size_t i = 0; i < stage.quests.size (); ++i) {
        const Quest &q = stage.quests[i];
        stage.cur[i] = distanceOf (m, q);
        // the tokens already on the place, up to what the quest needs, stay there
        if (m.get (q.place) > 0) {
          if (frozenNeed[q.place] == 0) frozenList.push_back (q.place);
          frozenNeed[q.place] = std::max<long long> (frozenNeed[q.place], std::min<long long> (q.need, m.get (q.place)));
        }
        total += stage.cur[i];
      }
      return total;
    }

    /** Would t drop a frozen place below the tokens a quest keeps there. */
    bool breaksFrozen (uint32_t t) const
    {
      if (frozenList.empty ()) return false;
      const SparseArray<T> &eff = net.effect (t);
      for (size_t i = 0; i < eff.size (); ++i) {
        size_t p = eff.keyAt (i);
        if (eff.valueAt (i) < 0 && frozenNeed[p] > 0 && current->get (p) + eff.valueAt (i) < frozenNeed[p]) return true;
      }
      return false;
    }

    /** Whether t moves a token of a component holding the quest's place (the tables of Components say so). */
    bool touches (uint32_t t, const Quest &q) const
    {
      for (const auto &cl : comps.consumedBy (t))
        if (comps.component (cl.first).indexOf (q.place) >= 0) return true;
      for (const auto &cl : comps.producedBy (t))
        if (comps.component (cl.first).indexOf (q.place) >= 0) return true;
      return false;
    }

    /** The stage's distance after firing t, recomputed for the quests t touches. */
    uint64_t score (uint32_t t, const Stage &stage, uint64_t total) const
    {
      const SparseArray<T> &eff = net.effect (t);
      uint64_t s = total;
      for (size_t i = 0; i < stage.quests.size (); ++i) {
        if (stage.cur[i] == 0) continue;
        const Quest &q = stage.quests[i];
        if (q.comp == NONE) {
          if (eff.get (q.place) > 0 && s > 0) s -= 1; // unguided: a transition feeding the place is the only clue
          continue;
        }
        if (!touches (t, q)) continue;
        s = s - stage.cur[i] + distanceOf (*current, q, &eff);
      }
      return s;
    }

    /** Where each process stands, for the alignment costs of a stage choice. */
    void locate (const Marking<T> &m)
    {
      standing.assign (comps.size (), {});
      for (uint32_t c = 0; c < comps.size (); ++c) {
        const auto &k = comps.component (c);
        for (uint32_t j = 0; j < k.size (); ++j)
          if (m.get (k.places[j]) > 0) standing[c].push_back (j);
      }
    }

    /**
     * The cost of aligning every process on the pre-places of t by their own
     * moves; a pre-place already holding its tokens costs nothing whatever its
     * component can do, UNREACHED when a process cannot get to its place alone.
     */
    uint64_t alignCost (uint32_t t) const
    {
      const SparseArray<T> &pre = net.pre (t);
      uint64_t align = 0;
      for (size_t i = 0; i < pre.size (); ++i) {
        size_t p = pre.keyAt (i);
        if (current->get (p) >= pre.valueAt (i)) continue;
        const auto &cs = comps.componentsOf (p);
        if (cs.empty ()) return UNREACHED;
        const std::vector<uint32_t> &dist = comps.questDistances (cs[0], p);
        uint64_t d = UNREACHED;
        for (uint32_t j : standing[cs[0]])
          if (dist[j] < d) d = dist[j];
        if (d >= Components<T>::BARRIER_OFFSET) return UNREACHED;
        align += d;
      }
      return align;
    }

    /** The stage's distance after t, from the located marking. */
    uint64_t stageAfter (const Stage &stage, uint32_t t) const
    {
      const SparseArray<T> &eff = net.effect (t);
      uint64_t s = 0;
      for (size_t i = 0; i < stage.quests.size (); ++i) {
        const Quest &q = stage.quests[i];
        if (q.comp == NONE || !touches (t, q)) s += stage.cur[i];
        else s += distanceOf (*current, q, &eff);
      }
      return s;
    }

    /**
     * Push a stage for the top: the barrier every process can reach alone,
     * not tabu at this level, whose outcome brings the top's quests closest
     * (the cheapest to align among equals). False when there is none.
     */
    bool push (const Marking<T> &m)
    {
      Stage &top = stack.back ();
      locate (m);
      uint32_t best = NONE;
      uint64_t bestAfter = std::numeric_limits<uint64_t>::max ();
      uint64_t bestAlign = std::numeric_limits<uint64_t>::max ();
      for (uint32_t t = 0; t < net.transitionCount (); ++t) {
        if (comps.syncDegreeOf (t) <= 1) continue;
        if (std::find (top.tabu.begin (), top.tabu.end (), t) != top.tabu.end ()) continue;
        uint64_t align = alignCost (t);
        if (align == UNREACHED) continue;
        uint64_t after = stageAfter (top, t);
        if (std::tie (after, align) < std::tie (bestAfter, bestAlign)) {
          bestAfter = after;
          bestAlign = align;
          best = t;
        }
      }
      if (best == NONE) return false;
      ++replans;
      if (debugSteps > 0)
        std::cerr << "stage " << stack.size () << ": " << (top.barrier == NONE ? std::string ("goal") : "t" + std::to_string (top.barrier))
            << " at " << top.best << " stalled, pushing t" << best << " after " << bestAfter << " align " << bestAlign
            << " tabu " << top.tabu.size () << std::endl;
      Stage sub { best, {}, {} };
      sub.parentBefore = 0;
      for (uint64_t d : top.cur) sub.parentBefore += d;
      const SparseArray<T> &pre = net.pre (best);
      for (size_t i = 0; i < pre.size (); ++i) addQuest (pre.keyAt (i), pre.valueAt (i), sub.quests);
      stack.push_back (std::move (sub));
      maxDepth = std::max<uint64_t> (maxDepth, stack.size ());
      return true;
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_COMPONENTSTRATEGY_H_ */
