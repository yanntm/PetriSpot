/*
 * QuestSweep.h
 *
 * The sweep as a sequence of quests: with no focus given, the strategy picks
 * an open target itself, the one whose places are nearest by the summed
 * component distance from where the processes stand (each thread takes the
 * target of its own rank, so threads chase different ones), and runs a
 * ComponentStrategy toward it; when the target is claimed, the quest stalls
 * or the run restarts, it picks the next. Every target satisfied on the way
 * is claimed by the walker, so a quest that reaches a shared pre-set claims
 * every transition on it. See algorithm.md.
 */
#ifndef PETRI_WALK_QUESTSWEEP_H_
#define PETRI_WALK_QUESTSWEEP_H_

#include <algorithm>
#include <cstdint>
#include <limits>
#include <memory>
#include <random>
#include <vector>

#include "expr/Expression.h"
#include "walk/ComponentStrategy.h"
#include "walk/Components.h"
#include "walk/RareStrategy.h"
#include "walk/Strategy.h"
#include "walk/TargetSet.h"

namespace petri::walk
{

template<typename T>
  class QuestSweep : public Strategy<T>
  {
    static constexpr uint32_t NONE = std::numeric_limits<uint32_t>::max ();

    const WalkNet<T> &net;
    const Components<T> &comps;
    const TargetSet<T> &targets;
    unsigned rank;              // this thread's rank among the threads sharing the sweep
    unsigned epsilonPercent;
    size_t sampleSize;
    uint64_t stallLimit;
    std::unique_ptr<ComponentStrategy<T>> quest;
    RareStrategy<T> filler;     // the move when no target can be quested
    uint32_t current = NONE;
    std::vector<std::vector<uint32_t>> standing;
    std::vector<std::pair<uint64_t, uint32_t>> ranked;

  public:
    /** Instrumentation. */
    uint64_t retargets = 0;
    uint64_t claimedOwn = 0;    // quests whose own target was claimed while they ran

    QuestSweep (const WalkNet<T> &n, const Components<T> &c, const TargetSet<T> &tgs, unsigned threadRank,
                unsigned epsilon, size_t sample, uint64_t stall)
        : net (n), comps (c), targets (tgs), rank (threadRank), epsilonPercent (epsilon), sampleSize (sample),
          stallLimit (stall), filler (sample)
    {
    }

    void onReset () override
    {
      // a restart may land near another target: choose again
      current = NONE;
      quest.reset ();
    }

    bool bestOfRun (SparseArray<T> &m, uint64_t &h) const override
    {
      return quest ? quest->bestOfRun (m, h) : false;
    }

    uint32_t choose (WalkContext<T> &ctx) override
    {
      if (current != NONE && targets.isSolved (current)) {
        ++claimedOwn;
        current = NONE;
      }
      if (current == NONE && !retarget (ctx.marking, ctx.rng)) return filler.choose (ctx);
      uint32_t t = quest->choose (ctx);
      if (t == RESTART) current = NONE; // the quest gave up: the walker restarts, the next choose picks again
      return t;
    }

  private:
    /** Where each process stands. */
    void locate (const Marking<T> &m)
    {
      standing.assign (comps.size (), {});
      for (uint32_t c = 0; c < comps.size (); ++c) {
        const auto &k = comps.component (c);
        for (uint32_t j = 0; j < k.size (); ++j)
          if (m.get (k.places[j]) > 0) standing[c].push_back (j);
      }
    }

    /** Summed component distance of the goal's atoms from the located marking; UNREACHED when an atom is not `place >= k`. */
    uint64_t distance (const petri::expr::Expression &e, const Marking<T> &m) const
    {
      using petri::expr::Expression;
      switch (e.kind) {
      case Expression::Kind::And: {
        uint64_t s = 0;
        for (const auto &c : e.children) {
          uint64_t d = distance (c, m);
          if (d == UNREACHED) return UNREACHED;
          s += d;
        }
        return s;
      }
      case Expression::Kind::Atom: {
        const auto &a = e.atom;
        if (a.terms.size () != 1 || a.terms[0].second <= 0 || a.op != petri::expr::Cmp::GE) return UNREACHED;
        size_t p = a.terms[0].first;
        long long need = (a.constant + a.terms[0].second - 1) / a.terms[0].second;
        if (m.get (p) >= need) return 0;
        const auto &cs = comps.componentsOf (p);
        if (cs.empty ()) return 1;
        const std::vector<uint32_t> &dist = comps.questDistances (cs[0], p);
        uint64_t d = UNREACHED;
        for (uint32_t j : standing[cs[0]])
          if (dist[j] < d) d = dist[j];
        return d == UNREACHED ? 1000000 : d + 1;
      }
      default:
        return UNREACHED;
      }
    }

    /** Open targets to rank at a retarget: all of them below this, a random sample above. */
    static constexpr size_t RANK_SAMPLE = 4096;

    /** Rank the open targets by distance and take the one of this thread's rank; false when none can be quested. */
    bool retarget (const Marking<T> &m, std::mt19937_64 &rng)
    {
      ++retargets;
      locate (m);
      ranked.clear ();
      size_t n = targets.size ();
      // a sample of the open targets when there are many: the nearest of a few thousand is near enough
      size_t stride = n > RANK_SAMPLE ? n / RANK_SAMPLE : 1;
      uint32_t first = stride > 1 ? static_cast<uint32_t> (rng () % stride) : 0;
      for (uint32_t id = first; id < n; id += static_cast<uint32_t> (stride)) {
        if (targets.isSolved (id)) continue;
        const Target<T> &tg = targets.target (id);
        if (tg.isDeadlock () || tg.isBound ()) continue;
        uint64_t d = distance (tg.expression (), m);
        if (d != UNREACHED) ranked.emplace_back (d, id);
      }
      if (ranked.empty () && stride > 1) {
        // the sample held no open target: look at them all
        for (uint32_t id = 0; id < n; ++id) {
          if (targets.isSolved (id)) continue;
          const Target<T> &tg = targets.target (id);
          if (tg.isDeadlock () || tg.isBound ()) continue;
          uint64_t d = distance (tg.expression (), m);
          if (d != UNREACHED) ranked.emplace_back (d, id);
        }
      }
      if (ranked.empty ()) return false;
      size_t pick = std::min<size_t> (rank, ranked.size () - 1);
      std::nth_element (ranked.begin (), ranked.begin () + pick, ranked.end ());
      current = ranked[pick].second;
      quest = std::make_unique<ComponentStrategy<T>> (net, comps, targets.target (current).expression (), epsilonPercent,
                                                       sampleSize, stallLimit);
      return true;
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_QUESTSWEEP_H_ */
