#pragma once
#include <algorithm>
#include <vector>
#include <cstddef>
namespace petri::reduction {
struct Graph {
  std::vector<std::vector<size_t>> next, previous;
  explicit Graph(size_t size) : next(size), previous(size) {}
  void edge(size_t from, size_t to) { next[from].push_back(to); }
  void finish() {
    for (size_t p = 0; p < next.size(); ++p) {
      auto& row = next[p];
      std::sort(row.begin(), row.end()); row.erase(std::unique(row.begin(), row.end()), row.end());
      for (size_t q : row) previous[q].push_back(p);
    }
  }
  template<class Stop> std::vector<std::vector<size_t>> components(Stop stop) const {
    std::vector<bool> seen(next.size(), false);
    std::vector<size_t> order;
    struct Frame { size_t node, edge; };
    std::vector<Frame> stack;
    for (size_t p = 0; p < next.size(); ++p) if (!seen[p]) {
      seen[p] = true; stack.push_back({p, 0});
      while (!stack.empty()) {
        if (stop()) return {};
        auto& f = stack.back();
        if (f.edge == next[f.node].size()) { order.push_back(f.node); stack.pop_back(); }
        else {
          size_t q = next[f.node][f.edge++];
          if (!seen[q]) { seen[q] = true; stack.push_back({q, 0}); }
        }
      }
    }
    std::fill(seen.begin(), seen.end(), false);
    std::vector<std::vector<size_t>> result;
    for (auto it = order.rbegin(); it != order.rend(); ++it) if (!seen[*it]) {
      result.emplace_back(); auto& component = result.back();
      component.push_back(*it); seen[*it] = true;
      for (size_t i = 0; i < component.size(); ++i) {
        if (stop()) return {};
        for (size_t q : previous[component[i]]) if (!seen[q]) { seen[q] = true; component.push_back(q); }
      }
    }
    return result;
  }
  template<class Stop> void prefix(std::vector<bool>& kept, Stop stop) const {
    std::vector<size_t> todo;
    for (size_t p = 0; p < kept.size(); ++p) if (kept[p]) todo.push_back(p);
    for (size_t i = 0; i < todo.size(); ++i) {
      if (stop()) return;
      for (size_t p : previous[todo[i]]) if (!kept[p]) { kept[p] = true; todo.push_back(p); }
    }
  }
};
}
