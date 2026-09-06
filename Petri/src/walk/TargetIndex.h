/*
 * TargetIndex.h
 *
 * A walker's index of the targets it checks: `own`, the ids it owns, and per
 * place the lists `up` and `down` of the own targets an increase or a
 * decrease of the place can make hold (from the mentions of TargetSet). A
 * firing touches a few places; the targets to re-evaluate are the up lists of
 * the places that grew and the down lists of those that shrank, merged once
 * each by an epoch stamp. Solved ids met while scanning leave the lists, so
 * they shrink over the run without a compaction pass; a walker whose own
 * targets are all solved takes over every open one. Deadlock targets are
 * checked apart and never indexed here. See algorithm.md.
 */
#ifndef PETRI_WALK_TARGETINDEX_H_
#define PETRI_WALK_TARGETINDEX_H_

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "walk/TargetSet.h"

namespace petri::walk
{

template<typename T>
  class TargetIndex
  {
    const TargetSet<T> &targets;
    std::vector<uint32_t> own;
    std::vector<std::vector<uint32_t>> up, down; // per place
    std::vector<uint32_t> stamp;                  // per target: epoch of its last candidacy
    std::vector<uint32_t> candidates;
    uint32_t epoch = 0;
    std::vector<typename TargetSet<T>::Mention> mentionBuf;

    /** Drop list[i] as solved; the caller does not advance. */
    static void dropAt (std::vector<uint32_t> &list, size_t i)
    {
      list[i] = list.back ();
      list.pop_back ();
    }

  public:
    /** Index the given ids, or every target when subset is null. */
    TargetIndex (const TargetSet<T> &tgs, const std::vector<uint32_t> *subset)
        : targets (tgs), up (tgs.places ()), down (tgs.places ()), stamp (tgs.size (), 0)
    {
      if (subset) {
        build (*subset);
      } else {
        std::vector<uint32_t> all (tgs.size ());
        for (uint32_t i = 0; i < all.size (); ++i) all[i] = i;
        build (all);
      }
    }

    void build (const std::vector<uint32_t> &ids)
    {
      own.clear ();
      for (auto &l : up) l.clear ();
      for (auto &l : down) l.clear ();
      for (uint32_t id : ids) {
        if (targets.target (id).isDeadlock () || targets.isSolved (id)) continue;
        own.push_back (id);
        targets.mentions (id, mentionBuf);
        for (const auto &m : mentionBuf) {
          if (m.direction >= 0) up[m.place].push_back (id);
          if (m.direction <= 0) down[m.place].push_back (id);
        }
      }
    }

    /** When every own target is solved, take over every open one. */
    void refill ()
    {
      if (own.empty () && targets.openCount () > targets.deadlocks ().size ()) build (targets.openTargets ());
    }

    /** f(id) on every own target still open; the solved ones met leave the list. */
    template<class F>
      void forEachOwn (F f)
      {
        for (size_t i = 0; i < own.size ();) {
          if (targets.isSolved (own[i])) { dropAt (own, i); continue; }
          f (own[i]);
          ++i;
        }
      }

    /** f(id) once on every own open target a touched place (place, grew) can have made hold. */
    template<class F>
      void forEachTouched (const std::vector<std::pair<size_t, bool>> &touched, F f)
      {
        ++epoch;
        candidates.clear ();
        for (const auto &pc : touched) {
          std::vector<uint32_t> &list = pc.second ? up[pc.first] : down[pc.first];
          for (size_t i = 0; i < list.size ();) {
            uint32_t id = list[i];
            if (targets.isSolved (id)) { dropAt (list, i); continue; }
            ++i;
            if (stamp[id] == epoch) continue;
            stamp[id] = epoch;
            candidates.push_back (id);
          }
        }
        for (uint32_t id : candidates) f (id);
      }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_TARGETINDEX_H_ */
