/*
 * NoveltyTracker.h
 *
 * The memory of one walker: how many times it fired each transition, which
 * ones it fired at least once, how long ago it last fired a new one, and the
 * marking reached by the last transition nobody had fired before (a rare
 * event, offered to the shared pool). Its counts are merged into the shared
 * Knowledge at every reset and every few thousand steps, over the
 * transitions touched since the last merge only. See algorithm.md.
 */
#ifndef PETRI_WALK_NOVELTYTRACKER_H_
#define PETRI_WALK_NOVELTYTRACKER_H_

#include <cstdint>
#include <vector>

#include "core/SparseArray.h"
#include "walk/Knowledge.h"
#include "walk/Marking.h"

namespace petri::walk
{

template<typename T>
  class NoveltyTracker
  {
    Knowledge *shared;               // may be null: the walker then only knows itself
    std::vector<uint32_t> fired;     // firings by this walker, per transition
    std::vector<uint32_t> pending;   // firings not yet merged into shared
    std::vector<uint32_t> dirty;     // transitions with a pending count
    std::vector<bool> firedOnce;
    uint64_t steps = 0;              // firings seen
    uint64_t lastNovel = 0;          // steps at the last first firing of a transition
    uint64_t distinct = 0;
    uint64_t rareEvents = 0;
    bool rarePending = false;        // the marking after a rare event awaits a snapshot
    bool hasRare = false;
    SparseArray<T> rareMarking;

  public:
    NoveltyTracker (size_t transitions, Knowledge *k)
        : shared (k), fired (transitions, 0), pending (transitions, 0), firedOnce (transitions, false)
    {
    }

    /** The walker fired t, times times; call before the marking moves on. */
    void onFired (uint32_t t, uint64_t times)
    {
      steps += times;
      bool novelHere = !firedOnce[t];
      if (novelHere) {
        firedOnce[t] = true;
        ++distinct;
        lastNovel = steps;
        // nobody fired it before, as far as the merged counts say
        if (fired[t] == 0 && (!shared || shared->fired (t) == 0)) {
          ++rareEvents;
          rarePending = true;
        }
      }
      fired[t] += static_cast<uint32_t> (times);
      if (pending[t] == 0) dirty.push_back (t);
      pending[t] += static_cast<uint32_t> (times);
    }

    /** After the firing: keep the marking reached by a rare event for the pool. */
    void observe (const Marking<T> &m)
    {
      if (!rarePending) return;
      rarePending = false;
      rareMarking = m.sparse ();
      hasRare = true;
    }

    /** Push the pending counts to the shared knowledge. */
    void merge ()
    {
      if (shared) {
        for (uint32_t t : dirty) shared->add (t, pending[t]);
      }
      for (uint32_t t : dirty) pending[t] = 0;
      dirty.clear ();
    }

    /** The marking of the last rare event, once; false when there is none. */
    bool takeRareMarking (SparseArray<T> &out)
    {
      if (!hasRare) return false;
      hasRare = false;
      out = rareMarking;
      return true;
    }

    /** Firings since a transition was last fired for the first time. */
    uint64_t stepsSinceNovelty () const
    {
      return steps - lastNovel;
    }
    /** Forget the wait, at a restart. */
    void resetNovelty ()
    {
      lastNovel = steps;
    }
    uint64_t distinctFired () const
    {
      return distinct;
    }
    uint64_t rareEventCount () const
    {
      return rareEvents;
    }
    const std::vector<uint32_t>& localCounts () const
    {
      return fired;
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_NOVELTYTRACKER_H_ */
