/*
 * TargetSet.h
 *
 * The open targets of a run, shared by the walker threads: names and
 * verdicts, atomic solved flags claimed by compare-and-swap, the deadlock
 * targets, and for each target the places it mentions with the direction of
 * change that can make it hold, from which every walker builds its own
 * index. See algorithm.md.
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
    size_t placeCount;
    std::vector<uint32_t> deadlockTargets;

  public:
    /** A place a target mentions, and whether an increase, a decrease or either can make the target hold. */
    struct Mention
    {
      size_t place;
      int direction; // +1 up, -1 down, 0 both
    };

  private:
    /**
     * The direction that can turn an atom true: for `sum >= c` an increase of
     * a place with a positive coefficient or a decrease of one with a negative
     * coefficient, the reverse for `<=`, either for `==` and `!=`; a negation
     * swaps the direction.
     */
    static void collectMentions (const petri::expr::Expression &e, bool negated, std::vector<Mention> &out)
    {
      using petri::expr::Cmp;
      if (e.kind == petri::expr::Expression::Kind::Atom) {
        for (const auto &t : e.atom.terms) {
          int dir;
          switch (e.atom.op) {
          case Cmp::GE: case Cmp::GT: dir = t.second > 0 ? 1 : -1; break;
          case Cmp::LE: case Cmp::LT: dir = t.second > 0 ? -1 : 1; break;
          default: dir = 0; break;
          }
          if (negated) dir = -dir;
          out.push_back ({ t.first, dir });
        }
        return;
      }
      bool neg = negated != (e.kind == petri::expr::Expression::Kind::Not);
      for (const auto &c : e.children) collectMentions (c, neg, out);
    }

    /** One entry per place, merged: two different directions give 0. */
    static void merge (std::vector<Mention> &ms)
    {
      std::sort (ms.begin (), ms.end (), [] (const Mention &a, const Mention &b) { return a.place < b.place; });
      size_t w = 0;
      for (size_t i = 0; i < ms.size (); ++i) {
        if (w > 0 && ms[w - 1].place == ms[i].place) {
          if (ms[w - 1].direction != ms[i].direction) ms[w - 1].direction = 0;
        } else {
          ms[w++] = ms[i];
        }
      }
      ms.resize (w);
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
          placeCount (placeCount)
    {
      open.store (targets.size ());
      for (uint32_t i = 0; i < targets.size (); ++i) {
        solved[i].store (false);
        best[i].store (std::numeric_limits<long long>::min ());
        if (targets[i].isDeadlock ()) deadlockTargets.push_back (i);
      }
    }

    size_t places () const
    {
      return placeCount;
    }

    /** The places target i mentions, each once, with the direction that can make i hold; a bound grows like a `>=`. */
    void mentions (uint32_t i, std::vector<Mention> &out) const
    {
      out.clear ();
      if (targets[i].isDeadlock ()) return;
      if (targets[i].isBound ()) {
        for (const auto &t : targets[i].boundForm ().terms) out.push_back ({ t.first, t.second > 0 ? 1 : -1 });
      } else {
        collectMentions (targets[i].expression (), false, out);
      }
      merge (out);
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
