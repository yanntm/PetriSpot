/*
 * RestartPolicy.h
 *
 * When a run ends and the walker restarts: a step budget, a wall clock
 * budget, a novelty stall, or any of them. The walker asks at every step
 * with a view of the run; the policy is built outside and injected, so a
 * schedule (Luby, adaptive) is one more class here. See algorithm.md.
 */
#ifndef PETRI_WALK_RESTARTPOLICY_H_
#define PETRI_WALK_RESTARTPOLICY_H_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace petri::walk
{

/** What a run looks like when the policy is asked. */
struct RunView
{
  uint64_t runSteps;          // steps since the run started
  uint64_t runMillis;         // wall clock since the run started (refreshed every 1024 steps)
  uint64_t stepsSinceNovelty; // steps since this walker last fired a transition for the first time
};

class RestartPolicy
{
public:
  virtual ~RestartPolicy () = default;
  virtual bool shouldRestart (const RunView &v) const = 0;
  /**
   * Whether the restart this policy asks for keeps the marking: the strategy
   * forgets (counters, tabu, current target), the walk goes on from where it
   * stands, a reachable state and, for a quest sweep, all of its progress.
   * False restores the initial marking or a pooled state.
   */
  virtual bool keepsMarking () const
  {
    return false;
  }
  /** What this policy is; for the report. */
  virtual std::string describe () const = 0;
};

/** The run ends after a number of steps; 0 never. */
class StepBudget : public RestartPolicy
{
  uint64_t steps;

public:
  explicit StepBudget (uint64_t n)
      : steps (n)
  {
  }
  bool shouldRestart (const RunView &v) const override
  {
    return steps != 0 && v.runSteps >= steps;
  }
  std::string describe () const override
  {
    return "steps " + std::to_string (steps);
  }
};

/** The run ends after a wall clock budget; 0 never. */
class WallTime : public RestartPolicy
{
  uint64_t millis;

public:
  explicit WallTime (uint64_t ms)
      : millis (ms)
  {
  }
  bool shouldRestart (const RunView &v) const override
  {
    return millis != 0 && v.runMillis >= millis;
  }
  /** The clock says the strategy has had its turn, not that the state is bad. */
  bool keepsMarking () const override
  {
    return true;
  }
  std::string describe () const override
  {
    return "wall " + std::to_string (millis) + " ms, in place";
  }
};

/** The run ends when no new transition was fired for a number of steps; 0 never. */
class NoveltyStall : public RestartPolicy
{
  uint64_t steps;

public:
  explicit NoveltyStall (uint64_t n)
      : steps (n)
  {
  }
  bool shouldRestart (const RunView &v) const override
  {
    return steps != 0 && v.stepsSinceNovelty >= steps;
  }
  std::string describe () const override
  {
    return "novelty stall " + std::to_string (steps);
  }
};

/** The run ends when any of the policies says so. */
class AnyOf : public RestartPolicy
{
  std::vector<std::unique_ptr<RestartPolicy>> policies;

public:
  void add (std::unique_ptr<RestartPolicy> p)
  {
    policies.push_back (std::move (p));
  }
  bool shouldRestart (const RunView &v) const override
  {
    for (const auto &p : policies)
      if (p->shouldRestart (v)) return true;
    return false;
  }
  /** As the first policy that asks for the restart says. */
  bool keepsMarking () const override
  {
    return false;
  }
  /** The verdict of the first policy asking for a restart on this view. */
  bool keepsMarking (const RunView &v) const
  {
    for (const auto &p : policies)
      if (p->shouldRestart (v)) return p->keepsMarking ();
    return false;
  }
  std::string describe () const override
  {
    std::string s;
    for (const auto &p : policies) s += (s.empty () ? "" : ", ") + p->describe ();
    return s;
  }
};

} // namespace petri::walk

#endif /* PETRI_WALK_RESTARTPOLICY_H_ */
