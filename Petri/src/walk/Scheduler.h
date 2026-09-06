/*
 * Scheduler.h
 *
 * Time sharing of many tasks over a few runner threads. Runnable tasks wait
 * in one queue ordered by virtual time; each runner pops the task with the
 * smallest clock, runs one slice (a number of steps, capped by the wall
 * clock), books the report, asks the decision hook whether the task goes on,
 * and pushes it back or finishes it. A task's clock advances by its running
 * time divided by its share, so shares translate into time. Tasks may be
 * added while the scheduler runs (a new task starts at the smallest clock in
 * the queue), and the spawn hook is asked for one whenever a runner would
 * otherwise idle. The run ends when no task is runnable and none can be
 * spawned, at the deadline, or on the stop flag; every task is then finished.
 * The policy (shares, parking, spawning) lives in the Coordinator; this file
 * is the queue and the runners. See WALK_PLAN.md sections 10.11 and 10.12.
 */
#ifndef PETRI_WALK_SCHEDULER_H_
#define PETRI_WALK_SCHEDULER_H_

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

#include "walk/Task.h"

namespace petri::walk
{

struct SchedulerSpec
{
  unsigned runners = 1;      // threads running the tasks
  uint64_t sliceSteps = 4096; // steps per slice
  uint64_t sliceMillis = 50;  // wall clock cap of a slice: a coarse step ends it
};

/** What the decision hook says of a task after a slice. */
enum class Decision
{
  Continue, Park
};

class Scheduler
{
public:
  /** Called under the scheduler's lock after each slice, with the task's index and the slice's report. */
  using OnSlice = std::function<Decision (size_t, const SliceReport&)>;
  /** Called under the lock when a runner would idle; null when nothing is worth spawning. */
  using Spawn = std::function<std::unique_ptr<Task> ()>;

private:
  std::vector<std::unique_ptr<Task>> tasks;
  SchedulerSpec spec;
  std::mutex mutex;
  std::condition_variable wake;
  std::vector<size_t> ready;   // indices of runnable tasks, a heap on vruntime (smallest first)
  size_t live = 0, running = 0; // tasks not finished; tasks in a slice right now
  bool halt = false;
  OnSlice onSlice;
  Spawn spawn;

  bool later (size_t a, size_t b) const
  {
    return tasks[a]->vruntime > tasks[b]->vruntime || (tasks[a]->vruntime == tasks[b]->vruntime && a > b);
  }
  void push (size_t i)
  {
    ready.push_back (i);
    std::push_heap (ready.begin (), ready.end (), [this] (size_t a, size_t b) { return later (a, b); });
  }
  size_t pop ()
  {
    std::pop_heap (ready.begin (), ready.end (), [this] (size_t a, size_t b) { return later (a, b); });
    size_t i = ready.back ();
    ready.pop_back ();
    return i;
  }

  /** Under the lock: add a task, its clock at the smallest in the queue so that it runs soon. */
  void admit (std::unique_ptr<Task> t)
  {
    t->vruntime = ready.empty () ? 0.0 : tasks[ready.front ()]->vruntime;
    tasks.push_back (std::move (t));
    ++live;
    push (tasks.size () - 1);
  }

  void runner (std::chrono::steady_clock::time_point deadline, const std::atomic<bool> *stop)
  {
    using clock = std::chrono::steady_clock;
    std::unique_lock<std::mutex> lock (mutex);
    for (;;) {
      // nothing runnable: spawn if the policy has something, else wait for a running task to come back
      while (!halt && ready.empty () && spawn) {
        std::unique_ptr<Task> t = spawn ();
        if (!t) break;
        admit (std::move (t));
      }
      while (!halt && ready.empty () && running > 0) wake.wait (lock);
      if (halt || ready.empty ()) break;
      if ((stop && stop->load (std::memory_order_relaxed)) || clock::now () >= deadline) {
        halt = true;
        wake.notify_all ();
        break;
      }
      size_t i = pop ();
      ++running;
      lock.unlock ();
      Slice slice { spec.sliceSteps, std::min (deadline, clock::now () + std::chrono::milliseconds (spec.sliceMillis)) };
      SliceReport r = tasks[i]->run (slice);
      lock.lock ();
      --running;
      tasks[i]->account (r);
      Decision d = onSlice ? onSlice (i, r) : Decision::Continue;
      if (r.finished || d == Decision::Park) {
        tasks[i]->finish ();
        --live;
      } else {
        push (i);
      }
      wake.notify_all ();
    }
  }

public:
  explicit Scheduler (SchedulerSpec s)
      : spec (s)
  {
  }

  void setHooks (OnSlice on, Spawn sp)
  {
    onSlice = std::move (on);
    spawn = std::move (sp);
  }

  /** Before run(): the initial tasks. */
  void add (std::unique_ptr<Task> t)
  {
    admit (std::move (t));
  }
  size_t size () const
  {
    return tasks.size ();
  }
  Task& task (size_t i)
  {
    return *tasks[i];
  }
  const Task& task (size_t i) const
  {
    return *tasks[i];
  }
  size_t liveCount () const
  {
    return live;
  }

  /** Run until nothing is runnable or spawnable, the deadline or the stop flag; then finish every task. */
  void run (std::chrono::steady_clock::time_point deadline, const std::atomic<bool> *stop)
  {
    running = 0;
    halt = false;
    unsigned n = std::max<unsigned> (1, spec.runners);
    if (n == 1) {
      runner (deadline, stop);
    } else {
      std::vector<std::thread> threads;
      for (unsigned k = 0; k < n; ++k) threads.emplace_back ([&] { runner (deadline, stop); });
      for (auto &th : threads) th.join ();
    }
    for (auto &t : tasks)
      if (!t->finished) t->finish ();
  }
};

} // namespace petri::walk

#endif /* PETRI_WALK_SCHEDULER_H_ */
