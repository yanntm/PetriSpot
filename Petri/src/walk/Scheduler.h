/*
 * Scheduler.h
 *
 * Time sharing of many tasks over a few runner threads. Runnable tasks wait
 * in one queue ordered by virtual time; each runner pops the task with the
 * smallest clock, runs one slice (a number of steps, capped by the wall
 * clock), books the report, and pushes the task back unless it finished. With
 * equal shares this is round robin; a task's clock advances by its running
 * time divided by its share, so shares translate into time. Shares follow
 * results: a slice's reward is its claims plus a small weight of the
 * transitions it fired for the first time; each strategy kind keeps its
 * reward and its running time both decayed exponentially with the running
 * time (a horizon of `tau` seconds), its score is the ratio, a reward rate
 * that no slice length can bias; the kinds' shares are their scores
 * normalised above a floor, split equally among the kind's tasks, so a kind
 * that stops paying keeps a trickle and comes back when it pays again. The
 * run ends when
 * every task has finished, the deadline has passed or the stop flag is set;
 * tasks still runnable are then left where they stand for `finish()`. See
 * WALK_PLAN.md section 10.11.
 */
#ifndef PETRI_WALK_SCHEDULER_H_
#define PETRI_WALK_SCHEDULER_H_

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
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
  bool adaptive = true;       // shares follow the reward rate per kind; false: equal shares throughout
  double shareFloor = 0.1;    // the part of the time shared equally whatever the scores (the rest follows them)
  double tau = 3.0;           // horizon in running seconds over which a kind's reward and time decay
  double noveltyWeight = 0.05; // a first firing counts this much of a claim
};

/** One strategy kind under the scheduler: its tasks, its score and its share of the time. */
struct KindShare
{
  std::string kind;
  unsigned tasks = 0;
  double reward = 0.0; // decayed reward of its tasks' slices
  double seconds = 0.0; // decayed running time of its tasks' slices
  double score = 0.0;  // reward / seconds: the reward rate over the horizon
  double share = 0.0;  // of the runners' time, the floor included
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
  std::vector<KindShare> kinds;
  std::vector<size_t> kindOf;  // per task

  /** Book a slice into its kind's score and redistribute the shares; under the lock. */
  void reward (size_t i, const SliceReport &r)
  {
    if (!spec.adaptive || r.micros == 0) return;
    KindShare &k = kinds[kindOf[i]];
    double seconds = static_cast<double> (r.micros) / 1e6;
    double decay = std::exp (-seconds / spec.tau);
    k.reward = k.reward * decay + static_cast<double> (r.claims) + spec.noveltyWeight * static_cast<double> (r.novelty);
    k.seconds = k.seconds * decay + seconds;
    k.score = k.seconds > 0 ? k.reward / k.seconds : 0.0;
    double total = 0.0;
    for (const KindShare &x : kinds) total += x.score;
    double equal = 1.0 / static_cast<double> (kinds.size ());
    for (KindShare &x : kinds)
      x.share = total > 0 ? spec.shareFloor * equal + (1.0 - spec.shareFloor) * x.score / total : equal;
    for (size_t t = 0; t < tasks.size (); ++t) {
      const KindShare &x = kinds[kindOf[t]];
      tasks[t]->share = x.share / static_cast<double> (x.tasks);
    }
  }

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
      reward (i, r);
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
    size_t k = 0;
    while (k < kinds.size () && kinds[k].kind != t->kind) ++k;
    if (k == kinds.size ()) kinds.push_back (KindShare { t->kind, 0, 0.0, 0.0, 0.0, 0.0 });
    ++kinds[k].tasks;
    kindOf.push_back (k);
    tasks.push_back (std::move (t));
  }

  /** The kinds with their final scores and shares. */
  const std::vector<KindShare>& shares () const
  {
    return kinds;
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
    for (KindShare &k : kinds) k.share = 1.0 / static_cast<double> (kinds.size ());
    for (size_t t = 0; t < tasks.size (); ++t) tasks[t]->share = kinds[kindOf[t]].share / kinds[kindOf[t]].tasks;
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
