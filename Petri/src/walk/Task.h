/*
 * Task.h
 *
 * An exploration task for the Scheduler: a resumable unit of work with a
 * label, a strategy kind, a share of the runners' time, a virtual clock and
 * its statistics. `run(slice)` advances it by at most `slice.steps` steps or
 * until `slice.deadline`, whichever comes first, and returns what the slice
 * produced; `finish()` closes it once the scheduler is done with it. A task
 * is owned by one runner at a time, so nothing in it locks. See
 * WALK_PLAN.md section 10.11.
 */
#ifndef PETRI_WALK_TASK_H_
#define PETRI_WALK_TASK_H_

#include <chrono>
#include <cstdint>
#include <string>

namespace petri::walk
{

struct Slice
{
  uint64_t steps;                                   // at most this many steps
  std::chrono::steady_clock::time_point deadline;   // and not past this instant (the cap on a coarse step)
};

struct SliceReport
{
  uint64_t steps = 0;      // steps taken in the slice
  uint64_t claims = 0;     // targets claimed in the slice
  uint64_t novelty = 0;    // transitions this task fired for the first time in the slice
  uint64_t micros = 0;     // running time of the slice
  bool capped = false;     // the deadline ended the slice before the steps did
  bool finished = false;   // the task has nothing left to do
};

class Task
{
public:
  std::string label;       // for the report, e.g. "sync:0:0" or "rare"
  std::string kind;        // the strategy name, what shares and summaries are grouped by
  double share = 1.0;      // weight in the scheduler: a task with twice the share gets twice the slices
  double vruntime = 0.0;   // virtual clock: advances by the slice's running time divided by the share
  // totals over the slices
  uint64_t steps = 0, micros = 0, claims = 0, novelty = 0, slices = 0, cappedSlices = 0, maxSliceMicros = 0;

  virtual ~Task () = default;
  virtual SliceReport run (const Slice &slice) = 0;
  virtual void finish () = 0;

  /** Book a slice's report into the totals and the virtual clock. */
  void account (const SliceReport &r)
  {
    steps += r.steps;
    micros += r.micros;
    claims += r.claims;
    novelty += r.novelty;
    ++slices;
    if (r.capped) ++cappedSlices;
    if (r.micros > maxSliceMicros) maxSliceMicros = r.micros;
    vruntime += static_cast<double> (r.micros) / (share > 0 ? share : 1e-9);
  }
};

} // namespace petri::walk

#endif /* PETRI_WALK_TASK_H_ */
