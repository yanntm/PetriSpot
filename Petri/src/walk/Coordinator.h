/*
 * Coordinator.h
 *
 * What decides after a slice. The coordinator owns the policy of the pool of
 * tasks: the shares of the strategy kinds (reward per running second, decayed
 * over a horizon), the parking of tasks that stopped progressing (their best
 * state goes to the shared pool, their record to the arm table), and the
 * spawning of new tasks from (state, tool) pairs chosen as a bandit: the
 * yield per running second of what ran from that pooled state with that
 * tool, optimistic when untried, so every tool runs from the initial marking
 * first and the estimates then concentrate the spawns. The catalogue is the
 * list of strategy specs; the factory the portfolio gives builds a task for
 * (tool, state). Bounded: the pool is a few dozen entries, the table pool
 * times catalogue. See WALK_PLAN.md section 10.12.
 */
#ifndef PETRI_WALK_COORDINATOR_H_
#define PETRI_WALK_COORDINATOR_H_

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "core/SparseArray.h"
#include "walk/Scheduler.h"
#include "walk/SharedPool.h"
#include "walk/Task.h"
#include "walk/WalkTask.h"

namespace petri::walk
{

struct CoordinatorSpec
{
  bool adaptive = true;        // shares follow the reward rate per kind; false: equal shares throughout
  double shareFloor = 0.1;     // the part of the time shared equally whatever the scores
  double tau = 3.0;            // horizon in running seconds over which a kind's reward and time decay
  double noveltyWeight = 0.05; // a first firing counts this much of a claim
  unsigned liveCap = 3;        // live tasks per runner
  unsigned grant = 8;          // slices a task may run without progress before it is parked
  uint64_t minRunMicros = 250000; // a task is not parked before it ran this long: no churn of tiny slices
  double armTau = 2.0;         // horizon in running seconds over which an arm's yield and time decay
  bool spawning = true;        // park and spawn; false: the initial tasks run to the end (step 1 and 2 behaviour)
};

/** One strategy kind: its tasks, its score and its share of the time. */
struct KindShare
{
  std::string kind;
  unsigned tasks = 0;
  double reward = 0.0, seconds = 0.0; // decayed over the horizon
  double score = 0.0;                 // reward / seconds
  double share = 0.0;                 // of the runners' time, the floor included
};

template<typename T>
  class Coordinator
  {
  public:
    /** Build a task of the tool from a state (null: the initial marking); the portfolio's factory. */
    using Factory = std::function<std::unique_ptr<WalkTask<T>> (size_t tool, uint64_t stateId, const SparseArray<T> *start)>;
    /** Whether anything is left to do; spawning stops when it says no. */
    using More = std::function<bool ()>;

  private:
    struct Arm
    {
      double yield = 0.0, seconds = 0.0; // decayed over armTau of the arm's own running time
      unsigned runs = 0;
    };
    struct Record
    {
      unsigned countdown = 0; // slices left without progress before parking
    };

    CoordinatorSpec spec;
    unsigned runners;
    std::vector<std::string> catalogue; // kind names, one per tool
    Factory factory;
    More more;
    SharedPool<T> *pool;
    std::vector<KindShare> kinds;       // per tool
    std::map<std::pair<uint64_t, size_t>, Arm> arms; // (state id, tool)
    std::vector<Record> records;        // per task index
    Scheduler *scheduler = nullptr;
    size_t nextTool = 0;                // round robin of the initial spawns
    uint64_t spawned = 0, parked = 0;
    SparseArray<T> scratch;

    bool progressed (const SliceReport &r) const
    {
      return r.claims > 0 || r.novelty > 0 || r.heuristicDrop;
    }

    /** Shares from the kinds' scores, split among each kind's live tasks. */
    void redistribute ()
    {
      double total = 0.0;
      for (const KindShare &k : kinds) total += k.score;
      double equal = 1.0 / static_cast<double> (kinds.size ());
      for (KindShare &k : kinds)
        k.share = (spec.adaptive && total > 0) ? spec.shareFloor * equal + (1.0 - spec.shareFloor) * k.score / total : equal;
      for (size_t i = 0; i < scheduler->size (); ++i) {
        Task &t = scheduler->task (i);
        if (t.finished) continue;
        const KindShare &k = kinds[t.tool];
        t.share = k.share / static_cast<double> (std::max<unsigned> (1, k.tasks));
      }
    }

    void book (size_t i, const SliceReport &r)
    {
      Task &t = scheduler->task (i);
      KindShare &k = kinds[t.tool];
      double seconds = static_cast<double> (r.micros) / 1e6;
      double decay = std::exp (-seconds / spec.tau);
      k.reward = k.reward * decay + static_cast<double> (r.claims) + spec.noveltyWeight * static_cast<double> (r.novelty);
      k.seconds = k.seconds * decay + seconds;
      k.score = k.seconds > 0 ? k.reward / k.seconds : 0.0;
    }

    /** A task ends (parked or finished): its yield goes to its arm, its state to the pool. */
    void close (size_t i, bool parkIt)
    {
      Task &t = scheduler->task (i);
      Arm &a = arms[{ t.stateId, t.tool }];
      double seconds = static_cast<double> (t.micros) / 1e6;
      double decay = std::exp (-seconds / spec.armTau);
      a.yield = a.yield * decay + static_cast<double> (t.claims) + spec.noveltyWeight * static_cast<double> (t.novelty);
      a.seconds = a.seconds * decay + seconds;
      ++a.runs;
      --kinds[t.tool].tasks;
      if (parkIt) {
        ++parked;
        auto &wt = static_cast<WalkTask<T>&> (t);
        wt.markParked ();
        uint64_t h = 0;
        if (pool && wt.bestState (scratch, h)) pool->publish (scratch, h, std::vector<uint32_t> ());
      }
    }

    /** The (state, tool) to spawn: untried pairs first (tools in turn, from the initial marking then the pool), then the best estimate. */
    bool choose (uint64_t &stateId, size_t &tool, SparseArray<T> &start, bool &fromInitial)
    {
      std::vector<typename SharedPool<T>::Entry> states;
      if (pool) states = pool->snapshot ();
      // untried from the initial marking, tools in turn
      for (size_t k = 0; k < catalogue.size (); ++k) {
        size_t tl = (nextTool + k) % catalogue.size ();
        if (arms.find ({ 0, tl }) == arms.end ()) {
          nextTool = tl + 1;
          stateId = 0;
          tool = tl;
          fromInitial = true;
          return true;
        }
      }
      // untried from a pooled state: the best state first
      std::sort (states.begin (), states.end (), [] (const auto &a, const auto &b) { return a.heuristic < b.heuristic; });
      for (const auto &e : states)
        for (size_t tl = 0; tl < catalogue.size (); ++tl)
          if (arms.find ({ e.id, tl }) == arms.end ()) {
            stateId = e.id;
            tool = tl;
            start = e.marking;
            fromInitial = false;
            return true;
          }
      // all tried: the best yield per second, with an exploration bonus for the rarely run
      double n = 0.0;
      for (const auto &kv : arms) n += kv.second.runs;
      double best = -1.0;
      bool found = false;
      for (const auto &kv : arms) {
        const Arm &a = kv.second;
        double estimate = a.yield / (a.seconds + 0.1) + 0.1 * std::sqrt (std::log (n + 1.0) / (a.runs + 1.0));
        if (estimate <= best) continue;
        uint64_t sid = kv.first.first;
        const SparseArray<T> *m = nullptr;
        if (sid != 0) {
          for (const auto &e : states)
            if (e.id == sid) m = &e.marking;
          if (!m) continue; // evicted meanwhile
        }
        best = estimate;
        stateId = sid;
        tool = kv.first.second;
        fromInitial = sid == 0;
        if (m) start = *m;
        found = true;
      }
      return found;
    }

  public:
    Coordinator (CoordinatorSpec s, unsigned runnerCount, std::vector<std::string> tools, Factory f, More m, SharedPool<T> *p)
        : spec (s), runners (runnerCount), catalogue (std::move (tools)), factory (std::move (f)), more (std::move (m)), pool (p)
    {
      for (const std::string &name : catalogue) kinds.push_back (KindShare { name, 0, 0.0, 0.0, 0.0, 0.0 });
    }

    /** Hook into the scheduler and seed it: one task per tool from the initial marking, up to the live cap. */
    void install (Scheduler &sched)
    {
      scheduler = &sched;
      sched.setHooks ([this] (size_t i, const SliceReport &r) { return onSlice (i, r); },
                      spec.spawning ? Scheduler::Spawn ([this] { return spawnOne (); }) : Scheduler::Spawn ());
      unsigned cap = std::max<unsigned> (1, spec.liveCap * runners);
      unsigned initial = spec.spawning ? std::min<unsigned> (cap, static_cast<unsigned> (catalogue.size ()) * runners)
                                       : std::max<unsigned> (1, 2 * runners);
      for (unsigned k = 0; k < initial; ++k) {
        size_t tool = k % catalogue.size ();
        admit (factory (tool, 0, nullptr), tool, 0);
      }
    }

    Decision onSlice (size_t i, const SliceReport &r)
    {
      book (i, r);
      redistribute ();
      if (r.finished) {
        close (i, false);
        return Decision::Continue;
      }
      if (!spec.spawning) return Decision::Continue;
      Record &rec = records[i];
      if (progressed (r)) rec.countdown = spec.grant;
      else if (rec.countdown > 0) --rec.countdown;
      if (rec.countdown == 0 && scheduler->task (i).micros >= spec.minRunMicros) {
        close (i, true);
        return Decision::Park;
      }
      return Decision::Continue;
    }

    const std::vector<KindShare>& shares () const
    {
      return kinds;
    }
    uint64_t spawnedCount () const
    {
      return spawned;
    }
    uint64_t parkedCount () const
    {
      return parked;
    }
    size_t armCount () const
    {
      return arms.size ();
    }

  private:
    void admit (std::unique_ptr<WalkTask<T>> t, size_t tool, uint64_t stateId)
    {
      t->tool = tool;
      t->stateId = stateId;
      t->kind = catalogue[tool];
      ++kinds[tool].tasks;
      records.push_back (Record { spec.grant });
      ++spawned;
      if (!arms.count ({ stateId, tool })) arms[{ stateId, tool }] = Arm ();
      scheduler->add (std::move (t));
      redistribute ();
    }

    std::unique_ptr<Task> spawnOne ()
    {
      if (scheduler->liveCount () >= spec.liveCap * runners) return nullptr;
      if (more && !more ()) return nullptr;
      uint64_t stateId;
      size_t tool;
      bool fromInitial;
      if (!choose (stateId, tool, scratch, fromInitial)) return nullptr;
      std::unique_ptr<WalkTask<T>> t = factory (tool, stateId, fromInitial ? nullptr : &scratch);
      t->tool = tool;
      t->stateId = stateId;
      t->kind = catalogue[tool];
      ++kinds[tool].tasks;
      records.push_back (Record { spec.grant });
      ++spawned;
      if (!arms.count ({ stateId, tool })) arms[{ stateId, tool }] = Arm ();
      return t;
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_COORDINATOR_H_ */
