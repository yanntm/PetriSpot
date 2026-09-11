#pragma once

#include <chrono>
#include <stdexcept>
#include <string_view>

namespace petri::reduction {

enum class Goal { NONE, DEADLOCK, REACHABILITY, SI_LTL, LTL, LIVENESS,
                  STATESPACE, LI_LTL, SI_CTL };

inline Goal parseGoal(std::string_view name) {
  if (name == "NONE") return Goal::NONE;
  if (name == "DEADLOCK") return Goal::DEADLOCK;
  if (name == "REACHABILITY") return Goal::REACHABILITY;
  if (name == "SI_LTL") return Goal::SI_LTL;
  if (name == "LTL") return Goal::LTL;
  if (name == "LIVENESS") return Goal::LIVENESS;
  if (name == "STATESPACE") return Goal::STATESPACE;
  if (name == "LI_LTL") return Goal::LI_LTL;
  if (name == "SI_CTL") return Goal::SI_CTL;
  throw std::invalid_argument("Unknown reduction goal");
}

struct Configuration {
  Goal goal = Goal::NONE;
  bool agglomeration = true;
  bool relevance = true;
  size_t maxPasses = 1000;
  size_t maxComposedArcs = 100000;
  size_t maxNameBytes = 1024;
  size_t implicitDepth = 5;
  std::chrono::milliseconds timeLimit {15000};
};

inline bool permitsAgglomeration(const Configuration& c) {
  return c.agglomeration && (c.goal == Goal::REACHABILITY || c.goal == Goal::DEADLOCK);
}

} // namespace petri::reduction
