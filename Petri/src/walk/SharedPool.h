/*
 * SharedPool.h
 *
 * Bounded, thread-safe pool of promising markings shared by the walkers of
 * a portfolio. A walker publishes the best state of a run (by its heuristic)
 * when the pool has room or the state beats the worst entry; a restarting
 * walker may draw an entry instead of the initial marking. Entries that
 * repeatedly fail to lead to an improvement are evicted, so states doomed by
 * an earlier wrong choice do not accumulate.
 */
#ifndef PETRI_WALK_SHAREDPOOL_H_
#define PETRI_WALK_SHAREDPOOL_H_

#include <cstdint>
#include <mutex>
#include <random>
#include <vector>

#include "core/SparseArray.h"

namespace petri::walk
{

template<typename T>
  class SharedPool
  {
  public:
    struct Entry
    {
      SparseArray<T> marking;
      uint64_t heuristic = 0;
      std::vector<uint32_t> trace; // from the initial marking; empty when not recorded
      unsigned failures = 0;
      uint64_t id = 0;
    };

  private:
    std::mutex mutex;
    std::vector<Entry> entries;
    size_t capacity;
    unsigned drawPercent;   // chance a restart draws from the pool
    unsigned maxFailures;   // evict after this many unproductive uses
    uint64_t nextId = 1;
    uint64_t published = 0, drawn = 0, evicted = 0;

  public:
    SharedPool (size_t cap, unsigned drawPct, unsigned maxFail = 3)
        : capacity (cap), drawPercent (drawPct), maxFailures (maxFail)
    {
    }

    /** Offer a state; returns true if it entered the pool. */
    bool publish (const SparseArray<T> &m, uint64_t h, const std::vector<uint32_t> &trace)
    {
      std::lock_guard<std::mutex> lock (mutex);
      size_t worst = 0;
      for (size_t i = 0; i < entries.size (); ++i) {
        if (entries[i].marking == m) return false;
        if (entries[i].heuristic > entries[worst].heuristic) worst = i;
      }
      if (entries.size () < capacity) {
        entries.push_back (Entry { m, h, trace, 0, nextId++ });
      } else if (h < entries[worst].heuristic) {
        entries[worst] = Entry { m, h, trace, 0, nextId++ };
        ++evicted;
      } else {
        return false;
      }
      ++published;
      return true;
    }

    /**
     * Maybe draw an entry for a restart; returns false to restart from the
     * initial marking. On success fills out and returns the entry id.
     */
    bool draw (std::mt19937_64 &rng, Entry &out)
    {
      std::lock_guard<std::mutex> lock (mutex);
      if (entries.empty () || rng () % 100 >= drawPercent) return false;
      out = entries[rng () % entries.size ()];
      ++drawn;
      return true;
    }

    /** Report whether a run started from entry id improved on its heuristic. */
    void report (uint64_t id, bool improved)
    {
      std::lock_guard<std::mutex> lock (mutex);
      for (size_t i = 0; i < entries.size (); ++i) {
        if (entries[i].id != id) continue;
        if (improved) {
          entries[i].failures = 0;
        } else if (++entries[i].failures >= maxFailures) {
          entries[i] = entries.back ();
          entries.pop_back ();
          ++evicted;
        }
        return;
      }
    }

    size_t size ()
    {
      std::lock_guard<std::mutex> lock (mutex);
      return entries.size ();
    }
    uint64_t publishedCount () const { return published; }
    uint64_t drawnCount () const { return drawn; }
    uint64_t evictedCount () const { return evicted; }
  };

} // namespace petri::walk

#endif /* PETRI_WALK_SHAREDPOOL_H_ */
