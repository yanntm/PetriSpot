/*
 * RelaxedPlanStrategy.h
 *
 * FF-style choice: compute the relaxed plan from the current marking and
 * fire one of its helpful transitions (enabled, in the plan). Falls back to a
 * uniform random enabled transition when no helpful transition is enabled,
 * and takes an epsilon share of random moves. Asks for a restart when the
 * heuristic value has not improved for stallLimit steps (0 disables).
 */
#ifndef PETRI_WALK_RELAXEDPLANSTRATEGY_H_
#define PETRI_WALK_RELAXEDPLANSTRATEGY_H_

#include <cstdint>
#include <iostream>
#include <vector>

#include "walk/RelaxedPlan.h"
#include "walk/Strategy.h"

namespace petri::walk
{

template<typename T>
  class RelaxedPlanStrategy : public Strategy<T>
  {
    RelaxedPlan<T> plan;
    unsigned epsilonPercent;
    uint64_t stallLimit;
    uint64_t sinceImprovement = 0;
    uint64_t bestThisRun = petri::expr::INFINITE_DISTANCE;
    std::vector<uint32_t> candidates;

  public:
    /** Print the first debugSteps decisions on stderr. */
    uint64_t debugSteps = 0;

    /** Instrumentation. */
    uint64_t minHeuristic = petri::expr::INFINITE_DISTANCE;
    uint64_t runsReachingMin = 0;
    uint64_t stepsWithoutHelpful = 0;
    SparseArray<T> bestMarking;

    RelaxedPlanStrategy (const WalkNet<T> &net, const petri::expr::Expression &goal,
                         unsigned epsilon = 5, uint64_t stall = 0)
        : plan (net, goal), epsilonPercent (epsilon), stallLimit (stall)
    {
    }

    void onReset () override
    {
      sinceImprovement = 0;
      bestThisRun = petri::expr::INFINITE_DISTANCE;
    }

    uint32_t choose (WalkContext<T> &ctx) override
    {
      const size_t n = ctx.enabled.size ();
      if (stallLimit && sinceImprovement >= stallLimit) return RESTART;
      if (epsilonPercent > 0 && ctx.rng () % 100 < epsilonPercent) {
        ++sinceImprovement;
        return ctx.enabled.at (ctx.rng () % n);
      }
      plan.compute (ctx.marking);
      uint64_t h = plan.heuristic ();
      ++sinceImprovement;
      if (h < bestThisRun) {
        bestThisRun = h;
        sinceImprovement = 0;
        if (h < minHeuristic) {
          minHeuristic = h;
          runsReachingMin = 1;
          bestMarking = ctx.marking.sparse ();
        } else if (h == minHeuristic) {
          ++runsReachingMin;
        }
      }
      candidates.clear ();
      for (uint32_t t : plan.helpfulTransitions ()) {
        if (ctx.enabled.isEnabled (t)) candidates.push_back (t);
      }
      if (debugSteps > 0) {
        --debugSteps;
        const auto &tn = ctx.net.getNet ().getTnames ();
        std::cerr << "h=" << h << " enabled=" << n << " plan=" << plan.plan ().size ()
            << " helpful=" << plan.helpfulTransitions ().size () << " [";
        for (uint32_t t : plan.helpfulTransitions ()) std::cerr << " " << tn[t] << (ctx.enabled.isEnabled (t) ? "" : "(disabled)");
        std::cerr << " ] plan: [";
        for (uint32_t t : plan.plan ()) std::cerr << " " << tn[t];
        std::cerr << " ]" << std::endl;
      }
      if (candidates.empty ()) {
        ++stepsWithoutHelpful;
        return ctx.enabled.at (ctx.rng () % n);
      }
      return candidates[ctx.rng () % candidates.size ()];
    }

    uint64_t initialHeuristic (const Marking<T> &m)
    {
      plan.compute (m);
      return plan.heuristic ();
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_RELAXEDPLANSTRATEGY_H_ */
