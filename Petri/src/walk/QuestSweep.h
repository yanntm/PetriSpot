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
 * every transition on it. A bound target on one place is a staircase: the
 * quest is `place >= best + 1`, and while the value climbs the same target is
 * quested again a step higher. See algorithm.md.
 */
#ifndef PETRI_WALK_QUESTSWEEP_H_
#define PETRI_WALK_QUESTSWEEP_H_

#include <algorithm>
#include <cstdint>
#include <limits>
#include <iostream>
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
    long long step = 0;         // bound targets: the value the current quest aims at
    petri::expr::Expression stepGoal; // its goal, kept alive for the quest
    std::vector<uint32_t> unquestable; // targets whose quest could not start this run
    std::vector<std::vector<uint32_t>> standing;
    std::vector<std::pair<uint64_t, uint32_t>> ranked;

  public:
    /** Print the first debugSteps retargets on stderr. */
    uint64_t debugSteps = 0;

    /** Instrumentation. */
    uint64_t retargets = 0;
    uint64_t claimedOwn = 0;    // quests whose own target was claimed while they ran
    uint64_t stepsClimbed = 0;  // bound targets raised by a step
    uint64_t abandoned = 0;     // quests that could not start or gave up
    uint64_t hopelessQuests = 0;    // of which: a quest's tokens could not reach its place
    uint64_t unstageableQuests = 0; // of which: behind barriers with no barrier to stage

    QuestSweep (const WalkNet<T> &n, const Components<T> &c, const TargetSet<T> &tgs, unsigned threadRank,
                unsigned epsilon, size_t sample, uint64_t stall)
        : net (n), comps (c), targets (tgs), rank (threadRank), epsilonPercent (epsilon), sampleSize (sample),
          stallLimit (stall), filler (sample)
    {
    }

    void onReset () override
    {
      // a restart may land near another target: choose again; what could not be quested stays so
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
      if (current != NONE && targets.target (current).isBound () && known (current, ctx.marking) >= step) {
        // the bound climbed a step: aim at the next one from here
        ++stepsClimbed;
        if (!aimBound (current, ctx.marking)) {
          unquestable.push_back (current);
          current = NONE;
        }
      }
      if (current == NONE && !retarget (ctx.marking, ctx.rng)) return filler.choose (ctx);
      uint32_t t = quest->choose (ctx);
      if (t == RESTART) {
        ++abandoned;
        hopelessQuests += quest->hopeless;
        unstageableQuests += quest->unstageable;
        if (quest->hopeless > 0) {
          // a process stands where nothing can reach its place: the marking is what is wrong, restart the walk
          current = NONE;
          return RESTART;
        }
        // behind barriers with none to stage from here: leave the target for this run, move on rarity
        unquestable.push_back (current);
        current = NONE;
        return filler.choose (ctx);
      }
      return t;
    }

  private:
    bool isUnquestable (uint32_t id) const
    {
      return std::find (unquestable.begin (), unquestable.end (), id) != unquestable.end ();
    }

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

    /** The best value known of a single-place bound: the shared record or the marking itself, whichever is higher. */
    long long known (uint32_t id, const Marking<T> &m) const
    {
      long long best = targets.bestValue (id);
      if (best == std::numeric_limits<long long>::min ()) best = 0;
      return std::max<long long> (best, m.get (targets.target (id).boundForm ().terms[0].first));
    }

    /** The next step of a single-place bound target as a goal, `place >= known + 1`; false for a sum or a reached limit. */
    bool aimBound (uint32_t id, const Marking<T> &m)
    {
      const Target<T> &tg = targets.target (id);
      const auto &form = tg.boundForm ();
      if (form.terms.size () != 1 || form.terms[0].second <= 0) return false;
      step = known (id, m) + 1;
      if (tg.hasLimit () && step > tg.getLimit ()) return false;
      petri::expr::LinearAtom a;
      a.addTerm (form.terms[0].first, 1);
      a.constant = step;
      a.op = petri::expr::Cmp::GE;
      stepGoal = petri::expr::Expression::makeAtom (std::move (a));
      quest = std::make_unique<ComponentStrategy<T>> (net, comps, stepGoal, epsilonPercent, sampleSize, stallLimit);
      return quest->supported ();
    }

    /** The goal of a target as an expression: its own for a reachability target, the next step of a single-place bound. */
    bool goalOf (uint32_t id, petri::expr::Expression &out, const Marking<T> &m) const
    {
      const Target<T> &tg = targets.target (id);
      if (tg.isDeadlock ()) return false;
      if (!tg.isBound ()) {
        out = tg.expression ();
        return true;
      }
      const auto &form = tg.boundForm ();
      if (form.terms.size () != 1 || form.terms[0].second <= 0) return false;
      long long next = known (id, m) + 1;
      if (tg.hasLimit () && next > tg.getLimit ()) return false;
      petri::expr::LinearAtom a;
      a.addTerm (form.terms[0].first, 1);
      a.constant = next;
      a.op = petri::expr::Cmp::GE;
      out = petri::expr::Expression::makeAtom (std::move (a));
      return true;
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
        long long have = m.get (p);
        if (have >= need) return 0;
        const auto &cs = comps.componentsOf (p);
        if (cs.empty ()) return UNGUIDED; // nothing to steer by: after every guided target
        // the missing tokens come from the components holding the place: the nearest ones
        std::vector<uint64_t> ds;
        for (uint32_t c : cs) {
          const auto &k = comps.component (c);
          const std::vector<uint32_t> &dist = comps.questDistances (c, p);
          for (uint32_t j : standing[c]) {
            if (k.places[j] == p || dist[j] == UNREACHED) continue;
            for (long long u = 0, v = m.get (k.places[j]); u < v; ++u) ds.push_back (dist[j]);
          }
        }
        size_t missing = static_cast<size_t> (need - have);
        if (ds.size () < missing) return HOPELESS;
        std::nth_element (ds.begin (), ds.begin () + static_cast<long> (missing) - 1, ds.end ());
        uint64_t sum = 0;
        for (size_t i = 0; i < missing; ++i) sum += ds[i] + 1;
        return sum;
      }
      default:
        return UNREACHED;
      }
    }

    /** Open targets to rank at a retarget: all of them below this, a random sample above. */
    static constexpr size_t RANK_SAMPLE = 4096;
    /** The rank distance of a target on a place no component holds: below a barrier, above any walk. */
    static constexpr uint64_t UNGUIDED = 900;
    /** The distance of a target whose tokens cannot reach it from the current marking. */
    static constexpr uint64_t HOPELESS = 1000000;

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
      petri::expr::Expression goal;
      for (uint32_t id = first; id < n; id += static_cast<uint32_t> (stride)) {
        if (targets.isSolved (id) || isUnquestable (id) || !goalOf (id, goal, m)) continue;
        uint64_t d = distance (goal, m);
        // a target no token can reach from here is not a quest from here
        if (d != UNREACHED && d < HOPELESS) ranked.emplace_back (d, id);
      }
      if (ranked.empty () && stride > 1) {
        // the sample held no open target: look at them all
        for (uint32_t id = 0; id < n; ++id) {
          if (targets.isSolved (id) || isUnquestable (id) || !goalOf (id, goal, m)) continue;
          uint64_t d = distance (goal, m);
          if (d != UNREACHED && d < HOPELESS) ranked.emplace_back (d, id);
        }
      }
      if (ranked.empty () && !unquestable.empty ()) {
        // everything left was given up on once: give it all another chance
        unquestable.clear ();
        return retarget (m, rng);
      }
      if (ranked.empty ()) return false;
      size_t pick = std::min<size_t> (rank, ranked.size () - 1);
      std::nth_element (ranked.begin (), ranked.begin () + pick, ranked.end ());
      current = ranked[pick].second;
      if (debugSteps > 0) {
        --debugSteps;
        std::cerr << "retarget: " << ranked.size () << " candidates, took " << targets.name (current) << " at distance "
            << ranked[pick].first << (targets.target (current).isBound () ? " (bound, known " + std::to_string (known (current, m)) + ")" : "")
            << std::endl;
      }
      if (targets.target (current).isBound ()) {
        if (aimBound (current, m)) return true;
        unquestable.push_back (current);
        current = NONE;
        return false;
      }
      quest = std::make_unique<ComponentStrategy<T>> (net, comps, targets.target (current).expression (), epsilonPercent,
                                                       sampleSize, stallLimit);
      return true;
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_QUESTSWEEP_H_ */
