/*
 * Scheduler.h
 *
 * Time sharing of many tasks over a few runner threads. Runnable tasks wait
 * in one queue ordered by virtual time; each runner pops the task with the
 * smallest clock, runs one slice (a number of steps, capped by the wall
 * clock), books the report, and pushes the task back unless it finished. With
 * equal shares this is round robin; a task's clock advances by its running
 * time divided by its share, so shares translate into time. The run ends when
 * every task has finished, the deadline has passed or the stop flag is set;
 * tasks still runnable are then left where they stand for `finish()`. See
 * WALK_PLAN.md section 10.11.
 */
#ifndef PETRI_WALK_SCHEDULER_H_
#define PETRI_WALK_SCHEDULER_H_

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
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
  unsigned tasks = 0;        // tasks to create (a portfolio setting; 0: as many as runners)
  uint64_t sliceSteps = 4096; // steps per slice
  uint64_t sliceMillis = 50;  // wall clock cap of a slice: a coarse step ends it
};

class Scheduler
{
  std::vector<std::unique_ptr<Task>> tasks;
  SchedulerSpec spec;
  std::mutex mutex;
  std::condition_variable wake;
  std::vector<size_t> ready;   // indices of runnable tasks, a heap on vruntime (smallest first)
  size_t finished = 0, running = 0;
  bool halt = false;

  static bool later (const std::vector<std::unique_ptr<Task>> &ts, size_t a, size_t b)
  {
    return ts[a]->vruntime > ts[b]->vruntime || (ts[a]->vruntime == ts[b]->vruntime && a > b);
  }

  void runner (std::chrono::steady_clock::time_point deadline, const std::atomic<bool> *stop)
  {
    using clock = std::chrono::steady_clock;
    std::unique_lock<std::mutex> lock (mutex);
    for (;;) {
      while (!halt && ready.empty () && finished + running < tasks.size ()) wake.wait (lock);
      if (halt || ready.empty ()) break;
      if ((stop && stop->load (std::memory_order_relaxed)) || clock::now () >= deadline) {
        halt = true;
        wake.notify_all ();
        break;
      }
      std::pop_heap (ready.begin (), ready.end (), [&] (size_t a, size_t b) { return later (tasks, a, b); });
      size_t i = ready.back ();
      ready.pop_back ();
      ++running;
      lock.unlock ();
      Slice slice { spec.sliceSteps, std::min (deadline, clock::now () + std::chrono::milliseconds (spec.sliceMillis)) };
      SliceReport r = tasks[i]->run (slice);
      lock.lock ();
      --running;
      tasks[i]->account (r);
      if (r.finished) {
        ++finished;
      } else {
        ready.push_back (i);
        std::push_heap (ready.begin (), ready.end (), [&] (size_t a, size_t b) { return later (tasks, a, b); });
      }
      wake.notify_one ();
    }
  }

public:
  explicit Scheduler (SchedulerSpec s)
      : spec (s)
  {
  }

  void add (std::unique_ptr<Task> t)
  {
    tasks.push_back (std::move (t));
  }
  size_t size () const
  {
    return tasks.size ();
  }
  Task& task (size_t i)
  {
    return *tasks[i];
  }

  /** Run every task under time sharing until all finish, the deadline or the stop flag; then finish them all. */
  void run (std::chrono::steady_clock::time_point deadline, const std::atomic<bool> *stop)
  {
    ready.clear ();
    for (size_t i = 0; i < tasks.size (); ++i) ready.push_back (i);
    std::make_heap (ready.begin (), ready.end (), [&] (size_t a, size_t b) { return later (tasks, a, b); });
    finished = 0;
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
    for (auto &t : tasks) t->finish ();
  }
};

} // namespace petri::walk

#endif /* PETRI_WALK_SCHEDULER_H_ */
