/*
 * WalkTask.h
 *
 * A Task over one Walker and its strategy: what one thread of the portfolio
 * ran before time sharing, made resumable. The first slice begins the walk
 * (from the initial marking or from the pooled state it was spawned on),
 * every slice advances it, finish() closes it into a ThreadReport and frees
 * the walker, so a parked task costs its report and nothing else. Claims are
 * published through the callback the portfolio gives, under its lock.
 */
#ifndef PETRI_WALK_WALKTASK_H_
#define PETRI_WALK_WALKTASK_H_

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "core/SparseArray.h"
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
  uint64_t stateId = 0; // the pooled state it started from, 0 for the initial marking
  bool parked = false;  // ended by the coordinator for want of progress
};

template<typename T>
  class WalkTask : public Task
  {
  public:
    /** Strategy specific counters and the best heuristic value, read at finish. */
    using Notes = std::function<void (ThreadReport&)>;

  private:
    std::unique_ptr<Strategy<T>> strategy; // owned here so that the walker's reference stays valid
    NoveltyTracker<T> tracker;
    std::unique_ptr<Walker<T>> walker;    // freed at finish
    WalkBudget budget;
    SparseArray<T> start;                 // the marking to begin from; empty for the initial marking
    bool fromState = false;
    Notes notes;
    bool begun = false;
    ThreadReport report;

  public:
    WalkTask (const WalkNet<T> &net, TargetSet<T> &targets, uint32_t focus, std::unique_ptr<Strategy<T>> st,
              uint64_t seed, Knowledge *knowledge, const RestartPolicy *policy, SharedPool<T> *pool, bool saturate,
              const WalkBudget &b, Notes n, const SparseArray<T> *startFrom = nullptr)
        : strategy (std::move (st)), tracker (net.transitionCount (), knowledge),
          walker (std::make_unique<Walker<T>> (net, targets, focus, *strategy, seed, nullptr)), budget (b), notes (std::move (n))
    {
      walker->setTracker (&tracker);
      walker->setKnowledge (knowledge);
      walker->setRestartPolicy (policy);
      walker->setPool (pool);
      walker->setSaturate (saturate);
      if (startFrom) {
        start = *startFrom;
        fromState = true;
      }
    }

    template<class F>
      void setOnClaim (F f)
      {
        walker->setOnClaim (std::move (f));
      }

    SliceReport run (const Slice &slice) override
    {
      if (!begun) {
        walker->begin (budget, fromState ? &start : nullptr);
        begun = true;
      }
      SliceReport r;
      walker->runSlice (slice.steps, slice.deadline, &r);
      return r;
    }

    /** The best state of the current run by the strategy's heuristic, for the pool when the task is parked. */
    bool bestState (SparseArray<T> &m, uint64_t &h) const
    {
      return walker && walker->bestState (m, h);
    }

    void finish () override
    {
      if (finished) return;
      if (!begun) walker->begin (budget, fromState ? &start : nullptr); // never scheduled: an empty walk, a consistent report
      WalkResult res = walker->finish ();
      report.strategy = label;
      report.stats = res.stats;
      report.found = res.found;
      report.claims = res.claims;
      report.slices = slices;
      report.cappedSlices = cappedSlices;
      report.maxSliceMicros = maxSliceMicros;
      report.stateId = stateId;
      if (notes) notes (report);
      walker.reset ();   // the dense working state goes; the strategy stays for the notes already taken
      strategy.reset ();
      finished = true;
    }

    void markParked ()
    {
      report.parked = true;
    }

    const ThreadReport& result () const
    {
      return report;
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_WALKTASK_H_ */
