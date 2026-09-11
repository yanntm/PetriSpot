// Supplementary finite-state check from module bring-up. MCC end-to-end
// validation in mcc.py is the primary acceptance path; see README.md.
#include <iostream>
#include <queue>
#include <random>
#include <set>
#include "reduction/Reduce.h"

using Net = SparsePetriNet<long>;
using State = std::vector<long>;
struct Exploration {
  std::set<State> states;
  bool deadlock = false;
  size_t edges = 0;
  long maxTotal = 0, maxPlace = 0;
};

static Exploration explore(const Net& net) {
  Exploration out;
  std::queue<State> todo;
  todo.push(net.getMarks()); out.states.insert(todo.front());
  while (!todo.empty()) {
    State state = todo.front(); todo.pop();
    bool enabled = false;
    long total = 0;
    for (long mark : state) { total += mark; out.maxPlace = std::max(out.maxPlace, mark); }
    out.maxTotal = std::max(out.maxTotal, total);
    for (size_t t = 0; t < net.getTransitionCount(); ++t) {
      const auto& pre = net.getFlowPT().getColumn(t);
      const auto& post = net.getFlowTP().getColumn(t);
      bool can = true;
      for (size_t i = 0; i < pre.size(); ++i) if (state[pre.keyAt(i)] < pre.valueAt(i)) can = false;
      if (!can) continue;
      enabled = true; ++out.edges;
      State next = state;
      for (size_t i = 0; i < pre.size(); ++i) next[pre.keyAt(i)] -= pre.valueAt(i);
      for (size_t i = 0; i < post.size(); ++i) next[post.keyAt(i)] += post.valueAt(i);
      if (out.states.insert(next).second) todo.push(std::move(next));
      if (out.states.size() > 100000) throw std::runtime_error("State limit");
    }
    out.deadlock |= !enabled;
  }
  return out;
}

static std::set<State> project(const Exploration& e, const std::vector<size_t>& coordinates) {
  std::set<State> result;
  for (const auto& state : e.states) {
    State selected;
    for (size_t p : coordinates) selected.push_back(state[p]);
    result.insert(std::move(selected));
  }
  return result;
}

int main(int argc, char** argv) {
  try {
    if (argc != 2) throw std::invalid_argument("Expected one model seed");
    unsigned seed = static_cast<unsigned>(std::stoul(argv[1]));
    std::mt19937 rng(seed);
    Net net;
    for (int p = 0; p < 6; ++p) net.addPlace("p" + std::to_string(p), p == 0 ? 3 : 0);
    // Non-increasing token total bounds this small state space.
    for (int t = 0; t < 9; ++t) {
      net.addTransition("t" + std::to_string(t));
      int a = t < 3 ? t : static_cast<int>(rng() % 6);
      int b = t < 3 ? t + 1 : static_cast<int>(rng() % 6);
      long w = t < 3 ? 1 : static_cast<long>(1 + rng() % 2);
      net.addPreArc(a, t, w);
      if (rng() % 5 != 0 || t < 3) net.addPostArc(b, t, w);
      if (t >= 3 && rng() % 3 == 0) {
        int read = static_cast<int>(rng() % 6);
        if (read != a && read != b) { net.addPreArc(read, t, 1); net.addPostArc(read, t, 1); }
      }
    }
    if (seed % 7 == 0) net.addTransition("source-identity");
    const auto original = explore(net);
    using namespace petri::reduction;
    for (Goal goal : {Goal::REACHABILITY, Goal::DEADLOCK, Goal::LTL, Goal::SI_CTL,
                      Goal::STATESPACE, Goal::LIVENESS, Goal::NONE}) {
      Configuration c; c.goal = goal;
      std::vector<bool> support(6, false);
      if (goal != Goal::DEADLOCK) { support[seed % 6] = true; support[(seed / 6) % 6] = true; }
      auto reduced = reduce(net, c, support);
      const auto after = explore(reduced.net);
      std::vector<size_t> oldCoordinates, newCoordinates;
      for (size_t p = 0; p < support.size(); ++p) if (support[p]) {
        oldCoordinates.push_back(p); newCoordinates.push_back(reduced.placeMap[p]);
      }
      if (project(original, oldCoordinates) != project(after, newCoordinates))
        throw std::runtime_error("Projected reachability mismatch");
      if (goal != Goal::REACHABILITY && original.deadlock != after.deadlock)
        throw std::runtime_error("Deadlock mismatch");
      if (goal == Goal::STATESPACE && (original.states.size() != after.states.size()
          || original.edges != after.edges || original.maxPlace != after.maxPlace || original.maxTotal != after.maxTotal))
        throw std::runtime_error("Counting mismatch");
    }
    std::cout << "VALID seed=" << seed << " states=" << original.states.size() << '\n';
  } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
