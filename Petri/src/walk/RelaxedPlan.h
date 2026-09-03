/*
 * RelaxedPlan.h
 *
 * Delete-relaxation heuristic over a marking (planning-style h_add):
 * a Dijkstra pass ignoring token consumption gives, for every place, the
 * additive cost of producing a token in it from the current marking;
 * transitions cost 1 plus the sum of their pre-place costs. The goal
 * expression is evaluated on these costs (And sums, Or takes the cheapest
 * branch) and a relaxed plan is extracted backward from the supporting goal
 * places; its transitions whose pre-places are all marked now are the
 * "helpful" transitions.
 */
#ifndef PETRI_WALK_RELAXEDPLAN_H_
#define PETRI_WALK_RELAXEDPLAN_H_

#include <algorithm>
#include <cstdint>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

#include "expr/Distance.h"
#include "expr/Expression.h"
#include "walk/Marking.h"
#include "walk/WalkNet.h"

namespace petri::walk
{

template<typename T>
  class RelaxedPlan
  {
    using Expression = petri::expr::Expression;
    static constexpr uint32_t INF = std::numeric_limits<uint32_t>::max ();
    static constexpr uint32_t NONE = std::numeric_limits<uint32_t>::max ();

    const WalkNet<T> &net;
    const Expression &goal;
    std::vector<size_t> goalPlaces;

    // per place; a stamp equal to epoch means the entry is valid this round
    std::vector<uint32_t> cost;
    std::vector<uint32_t> achiever;
    std::vector<uint32_t> placeStamp;
    std::vector<char> goalPending;
    // per transition
    std::vector<uint32_t> preSize;
    std::vector<uint32_t> remaining;
    std::vector<uint32_t> transitionStamp;
    std::vector<char> inPlan;
    std::vector<uint32_t> planTransitions;
    uint32_t epoch = 0;

    std::vector<size_t> supporters;
    std::vector<uint32_t> helpful;
    uint64_t value = 0;

    void collectGoalPlaces (const Expression &e)
    {
      if (e.kind == Expression::Kind::Atom) {
        const auto &a = e.atom;
        if (a.terms.size () == 1 && a.terms[0].second > 0 && a.constant >= 1
            && (a.op == petri::expr::Cmp::GE || a.op == petri::expr::Cmp::EQ)) {
          size_t p = a.terms[0].first;
          for (size_t q : goalPlaces) if (q == p) return;
          goalPlaces.push_back (p);
        }
        return;
      }
      for (const auto &c : e.children) collectGoalPlaces (c);
    }

    /** Cost of p this round (INF when never reached). */
    uint32_t costAt (size_t p) const
    {
      return placeStamp[p] == epoch ? cost[p] : INF;
    }
    void touchPlace (size_t p)
    {
      if (placeStamp[p] != epoch) {
        placeStamp[p] = epoch;
        cost[p] = INF;
        achiever[p] = NONE;
      }
    }

    void clear ()
    {
      ++epoch;
      if (epoch == 0) { // wrapped: invalidate everything explicitly
        std::fill (placeStamp.begin (), placeStamp.end (), 0);
        std::fill (transitionStamp.begin (), transitionStamp.end (), 0);
        epoch = 1;
      }
      for (uint32_t t : planTransitions) inPlan[t] = 0;
      planTransitions.clear ();
      supporters.clear ();
      helpful.clear ();
    }

    /** Dijkstra from the marked places until every goal place has a cost. */
    void propagate (const Marking<T> &m)
    {
      using Item = std::pair<uint32_t, size_t>;
      std::priority_queue<Item, std::vector<Item>, std::greater<Item>> queue;
      const SparseArray<T> &sp = m.sparse ();
      for (size_t i = 0; i < sp.size (); ++i) {
        size_t p = sp.keyAt (i);
        touchPlace (p);
        cost[p] = 0;
        queue.push ({ 0, p });
      }
      size_t goalsLeft = 0;
      for (size_t g : goalPlaces) if (costAt (g) == INF) ++goalsLeft;
      std::vector<char> &pending = goalPending;
      for (size_t g : goalPlaces) pending[g] = costAt (g) == INF ? 1 : 0;
      const MatrixCol<T> &post = net.getNet ().getFlowTP ();
      while (!queue.empty () && goalsLeft > 0) {
        auto [c, p] = queue.top ();
        queue.pop ();
        if (c > cost[p]) continue;
        if (pending[p]) {
          pending[p] = 0;
          --goalsLeft;
        }
        for (const auto &cons : net.consumersOf (p)) {
          uint32_t t = cons.transition;
          if (transitionStamp[t] != epoch) {
            transitionStamp[t] = epoch;
            remaining[t] = preSize[t];
          }
          if (--remaining[t] != 0) continue;
          uint64_t tc = 1;
          const SparseArray<T> &pre = net.pre (t);
          for (size_t j = 0; j < pre.size (); ++j) tc += cost[pre.keyAt (j)];
          uint32_t tcost = tc >= INF ? INF - 1 : static_cast<uint32_t> (tc);
          const SparseArray<T> &out = post.getColumn (t);
          for (size_t j = 0; j < out.size (); ++j) {
            size_t q = out.keyAt (j);
            touchPlace (q);
            if (tcost < cost[q]) {
              cost[q] = tcost;
              achiever[q] = t;
              queue.push ({ tcost, q });
            }
          }
        }
      }
    }

    /** h_add of e on the current costs; records the supporting goal places. */
    uint64_t evaluate (const Expression &e, const Marking<T> &m, bool record)
    {
      switch (e.kind) {
      case Expression::Kind::True: return 0;
      case Expression::Kind::False: return petri::expr::INFINITE_DISTANCE;
      case Expression::Kind::Not: return e.eval (m) ? 0 : 1;
      case Expression::Kind::And: {
        uint64_t sum = 0;
        for (const auto &c : e.children) {
          uint64_t d = evaluate (c, m, record);
          if (d >= petri::expr::INFINITE_DISTANCE) return petri::expr::INFINITE_DISTANCE;
          sum += d;
        }
        return sum;
      }
      case Expression::Kind::Or: {
        uint64_t best = petri::expr::INFINITE_DISTANCE;
        const Expression *bestChild = nullptr;
        for (const auto &c : e.children) {
          uint64_t d = evaluate (c, m, false);
          if (d < best) {
            best = d;
            bestChild = &c;
          }
        }
        if (record && bestChild) evaluate (*bestChild, m, true);
        return best;
      }
      case Expression::Kind::Atom: {
        const auto &a = e.atom;
        uint64_t md = petri::expr::atomDistance (a.value (m), a.op, a.constant);
        if (md == 0) return 0;
        if (a.terms.size () == 1 && a.terms[0].second > 0 && a.constant >= 1
            && (a.op == petri::expr::Cmp::GE || a.op == petri::expr::Cmp::EQ)) {
          size_t p = a.terms[0].first;
          uint32_t c = costAt (p);
          if (c == INF) return petri::expr::INFINITE_DISTANCE;
          if (record) supporters.push_back (p);
          return c > md ? c : md;
        }
        return md;
      }
      }
      return petri::expr::INFINITE_DISTANCE;
    }

    /** Backward extraction from the supporters; fills planTransitions and helpful. */
    void extract ()
    {
      std::vector<size_t> stack (supporters.begin (), supporters.end ());
      while (!stack.empty ()) {
        size_t g = stack.back ();
        stack.pop_back ();
        uint32_t cg = costAt (g);
        if (cg == 0 || cg == INF) continue;
        uint32_t t = achiever[g];
        if (t == NONE || inPlan[t]) continue;
        inPlan[t] = 1;
        planTransitions.push_back (t);
        bool applicable = true;
        const SparseArray<T> &pre = net.pre (t);
        for (size_t j = 0; j < pre.size (); ++j) {
          size_t r = pre.keyAt (j);
          if (costAt (r) != 0) {
            applicable = false;
            stack.push_back (r);
          }
        }
        if (applicable) helpful.push_back (t);
      }
    }

  public:
    RelaxedPlan (const WalkNet<T> &n, const Expression &g)
        : net (n), goal (g), cost (n.placeCount (), INF), achiever (n.placeCount (), NONE),
          placeStamp (n.placeCount (), 0), goalPending (n.placeCount (), 0),
          preSize (n.transitionCount (), 0), remaining (n.transitionCount (), 0),
          transitionStamp (n.transitionCount (), 0), inPlan (n.transitionCount (), 0)
    {
      collectGoalPlaces (goal);
      for (size_t t = 0; t < n.transitionCount (); ++t) {
        preSize[t] = static_cast<uint32_t> (n.pre (t).size ());
      }
    }

    size_t goalPlaceCount () const
    {
      return goalPlaces.size ();
    }

    /** Recompute costs, the heuristic value and the helpful transitions for m. */
    void compute (const Marking<T> &m)
    {
      clear ();
      propagate (m);
      value = evaluate (goal, m, true);
      extract ();
    }

    uint64_t heuristic () const
    {
      return value;
    }
    const std::vector<uint32_t>& helpfulTransitions () const
    {
      return helpful;
    }
    const std::vector<uint32_t>& plan () const
    {
      return planTransitions;
    }
    uint32_t costOf (size_t place) const
    {
      return costAt (place);
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_RELAXEDPLAN_H_ */
