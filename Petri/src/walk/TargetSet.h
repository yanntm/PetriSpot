/*
 * TargetSet.h
 *
 * The open targets of a run, shared by the walker threads: names and
 * verdicts, atomic solved flags claimed by compare-and-swap, the
 * place-to-targets index used for incremental checks, and the deadlock
 * targets. See algorithm.md.
 */
#ifndef PETRI_WALK_TARGETSET_H_
#define PETRI_WALK_TARGETSET_H_

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "expr/Expression.h"
#include "walk/Target.h"

namespace petri::walk
{

/** Focus value of a walker that aims at no particular target. */
constexpr uint32_t NO_FOCUS = std::numeric_limits<uint32_t>::max ();

template<typename T>
  class TargetSet
  {
    std::vector<Target<T>> targets;
    std::vector<std::string> names;
    std::vector<std::string> verdicts;
    std::unique_ptr<std::atomic<bool>[]> solved;
    std::unique_ptr<std::atomic<long long>[]> best; // bound targets: largest value seen
    std::atomic<size_t> open { 0 };
    std::vector<std::vector<uint32_t>> byPlace; // targets whose goal mentions the place
    std::vector<uint32_t> deadlockTargets;

    static void collectPlaces (const petri::expr::Expression &e, std::vector<size_t> &out)
    {
      if (e.kind == petri::expr::Expression::Kind::Atom) {
        for (const auto &t : e.atom.terms) out.push_back (t.first);
      }
      for (const auto &c : e.children) collectPlaces (c, out);
    }

  public:
    /**
     * Build from the targets with their names and MCC verdict words; the set
     * is fixed afterwards (the flags are atomics).
     */
    TargetSet (size_t placeCount, std::vector<Target<T>> tgs, std::vector<std::string> nms,
               std::vector<std::string> vds)
        : targets (std::move (tgs)), names (std::move (nms)), verdicts (std::move (vds)),
          solved (new std::atomic<bool>[targets.size ()]), best (new std::atomic<long long>[targets.size ()]),
          byPlace (placeCount)
    {
      open.store (targets.size ());
      std::vector<size_t> places;
      for (uint32_t i = 0; i < targets.size (); ++i) {
        solved[i].store (false);
        best[i].store (std::numeric_limits<long long>::min ());
        if (targets[i].isDeadlock ()) {
          deadlockTargets.push_back (i);
          continue;
        }
        places.clear ();
        if (targets[i].isBound ()) {
          for (const auto &t : targets[i].boundForm ().terms) places.push_back (t.first);
        } else {
          collectPlaces (targets[i].expression (), places);
        }
        std::sort (places.begin (), places.end ());
        places.erase (std::unique (places.begin (), places.end ()), places.end ());
        for (size_t p : places) byPlace[p].push_back (i);
      }
    }

    size_t size () const
    {
      return targets.size ();
    }
    size_t openCount () const
    {
      return open.load (std::memory_order_relaxed);
    }
    bool isSolved (uint32_t i) const
    {
      return solved[i].load (std::memory_order_relaxed);
    }
    /** Mark i solved; true iff this call did it (exactly one caller wins). */
    bool claim (uint32_t i)
    {
      bool expected = false;
      if (!solved[i].compare_exchange_strong (expected, true)) return false;
      open.fetch_sub (1, std::memory_order_relaxed);
      return true;
    }
    const Target<T>& target (uint32_t i) const
    {
      return targets[i];
    }
    /** Bound targets: raise the shared maximum to v if larger. */
    void recordValue (uint32_t i, long long v)
    {
      long long cur = best[i].load (std::memory_order_relaxed);
      while (v > cur && !best[i].compare_exchange_weak (cur, v, std::memory_order_relaxed)) {
      }
    }
    /** Bound targets: the largest value seen so far (LLONG_MIN before any). */
    long long bestValue (uint32_t i) const
    {
      return best[i].load (std::memory_order_relaxed);
    }
    const std::string& name (uint32_t i) const
    {
      return names[i];
    }
    const std::string& verdict (uint32_t i) const
    {
      return verdicts[i];
    }
    const std::vector<uint32_t>& targetsOf (size_t place) const
    {
      return byPlace[place];
    }
    const std::vector<uint32_t>& deadlocks () const
    {
      return deadlockTargets;
    }
    /** Indices of the targets not solved yet, in order. */
    std::vector<uint32_t> openTargets () const
    {
      std::vector<uint32_t> out;
      for (uint32_t i = 0; i < targets.size (); ++i)
        if (!isSolved (i)) out.push_back (i);
      return out;
    }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_TARGETSET_H_ */
