/*
 * Knowledge.h
 *
 * What the walker threads know together: how many times each transition has
 * been fired by any of them. Bounded (one counter per transition), updated by
 * relaxed atomics, read without locking; a stale value costs nothing but a
 * slightly worse choice. The rarity of a transition, 1 / (1 + fired), is what
 * the rare-event choice and the pool feed on. See algorithm.md.
 */
#ifndef PETRI_WALK_KNOWLEDGE_H_
#define PETRI_WALK_KNOWLEDGE_H_

#include <atomic>
#include <cstdint>
#include <memory>

namespace petri::walk
{

class Knowledge
{
  std::unique_ptr<std::atomic<uint32_t>[]> firedAll;
  size_t n;

public:
  explicit Knowledge (size_t transitions)
      : firedAll (new std::atomic<uint32_t>[transitions]), n (transitions)
  {
    for (size_t t = 0; t < n; ++t) firedAll[t].store (0, std::memory_order_relaxed);
  }

  size_t size () const
  {
    return n;
  }
  /** Firings of t by every thread, as merged so far. */
  uint32_t fired (uint32_t t) const
  {
    return firedAll[t].load (std::memory_order_relaxed);
  }
  /** Merge k firings of t from one thread. */
  void add (uint32_t t, uint32_t k)
  {
    firedAll[t].fetch_add (k, std::memory_order_relaxed);
  }
  /** Transitions fired at least once by anyone. */
  size_t distinctFired () const
  {
    size_t d = 0;
    for (size_t t = 0; t < n; ++t) d += fired (static_cast<uint32_t> (t)) > 0;
    return d;
  }
};

} // namespace petri::walk

#endif /* PETRI_WALK_KNOWLEDGE_H_ */
