/*
 * Verdict.h
 *
 * The three-valued verdict of a configuration, the budgets of a round and
 * the counters of a check.
 */
#ifndef PETRI_CTL_VERDICT_H_
#define PETRI_CTL_VERDICT_H_

#include <chrono>
#include <cstdint>
#include <ostream>

namespace petri::ctl
{

enum class Verdict : int8_t
{
  False = 0, True = 1, Unknown = -1
};

inline Verdict negate (Verdict v)
{
  return v == Verdict::True ? Verdict::False : v == Verdict::False ? Verdict::True : Verdict::Unknown;
}

inline const char* to_string (Verdict v)
{
  return v == Verdict::True ? "TRUE" : v == Verdict::False ? "FALSE" : "UNKNOWN";
}

struct Budget
{
  uint64_t huntSteps = 1000;     // steps of one hunt at the root of the search, over its restarts
  uint64_t regionStates = 1000;  // states a region at the root may hold
  uint64_t nestedSteps = 1000;   // the same for a hunt opened from a state met along a hunt or in a region
  uint64_t nestedStates = 1000;  // the same for a nested region: the probe of one state
  uint64_t runLength = 1000;     // steps of one run of a hunt before a restart
  unsigned epsilon = 10;         // percentage of uniform moves in a hunt
  bool saturate = true;          // every other run of a hunt repeats its chosen transition while it stays enabled
  size_t sample = 8;             // successors scored per greedy step (0: all)
  unsigned level = 0;            // the round: memo entries left unknown at a lower level are retried
  std::chrono::steady_clock::time_point deadline = std::chrono::steady_clock::time_point::max ();
};

struct Stats
{
  uint64_t solves = 0;        // solve calls on temporal nodes
  uint64_t memoHits = 0;
  uint64_t configurations = 0; // distinct (marking, node) pairs decided or tried
  uint64_t hunts = 0;
  uint64_t huntRuns = 0;
  uint64_t huntSteps = 0;
  uint64_t regions = 0;
  uint64_t regionStates = 0;
  uint64_t propagated = 0;     // verdicts remembered along witness paths and closed regions

  void print (std::ostream &os) const
  {
    os << "configurations " << configurations << ", memo hits " << memoHits << ", hunts " << hunts << " (" << huntRuns
        << " runs, " << huntSteps << " steps), regions " << regions << " (" << regionStates << " states), propagated "
        << propagated;
  }
};

} // namespace petri::ctl

#endif /* PETRI_CTL_VERDICT_H_ */
