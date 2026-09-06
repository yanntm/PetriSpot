/*
 * WalkTask.h
 *
 * A Task over one Walker and its strategy: what one thread of the portfolio
 * ran before time sharing, made resumable. The first slice begins the walk,
 * every slice advances it, finish() closes it into a ThreadReport. Claims
 * are published through the callback the portfolio gives, under its lock.
 */
#ifndef PETRI_WALK_WALKTASK_H_
#define PETRI_WALK_WALKTASK_H_

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "walk/NoveltyTracker.h"
#include "walk/Task.h"
#include "walk/Walker.h"

namespace petri::walk
{

/** What one task reports at the end: its strategy, counters and claims. */
struct ThreadReport
{
  std::string strategy;
  std::string notes;
  WalkStats stats;
  uint64_t minHeuristic = 0;
  bool found = false;   // claimed the focus
  uint64_t claims = 0;
  uint64_t slices = 0, cappedSlices = 0, maxSliceMicros = 0; // from the scheduler
};

template<typename T>
  class WalkTask : public Task
  {
  public:
    /** Strategy specific counters and the best heuristic value, read at finish. */
    using Notes = std::function<void (ThreadReport&)>;

  private:
    std::unique_ptr<Strategy<T>> strategy; // owned here so that the walker's reference stays valid
    std::vector<uint32_t> subset;         // the targets this task checks, empty for all
    NoveltyTracker<T> tracker;
    Walker<T> walker;
    WalkBudget budget;
    Notes notes;
    bool begun = false;
    ThreadReport report;

  public:
    WalkTask (const WalkNet<T> &net, TargetSet<T> &targets, uint32_t focus, std::unique_ptr<Strategy<T>> st,
              uint64_t seed, std::vector<uint32_t> owned, Knowledge *knowledge, const RestartPolicy *policy,
              SharedPool<T> *pool, bool saturate, const WalkBudget &b, Notes n)
        : strategy (std::move (st)), subset (std::move (owned)), tracker (net.transitionCount (), knowledge),
          walker (net, targets, focus, *strategy, seed, subset.empty () ? nullptr : &subset), budget (b), notes (std::move (n))
    {
      walker.setTracker (&tracker);
      walker.setKnowledge (knowledge);
      walker.setRestartPolicy (policy);
      walker.setPool (pool);
      walker.setSaturate (saturate);
    }

    template<class F>
      void setOnClaim (F f)
      {
        walker.setOnClaim (std::move (f));
      }

    SliceReport run (const Slice &slice) override
    {
      if (!begun) {
        walker.begin (budget);
        begun = true;
      }
      SliceReport r;
      walker.runSlice (slice.steps, slice.deadline, &r);
      return r;
    }

    void finish () override
    {
      if (!begun) walker.begin (budget); // never scheduled: an empty walk, so that the report is consistent
      WalkResult res = walker.finish ();
      report.strategy = label;
      report.stats = res.stats;
      report.found = res.found;
      report.claims = res.claims;
      report.slices = slices;
      report.cappedSlices = cappedSlices;
      report.maxSliceMicros = maxSliceMicros;
      if (notes) notes (report);
    }

    const ThreadReport& result () const
    {
      return report;
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_WALKTASK_H_ */
