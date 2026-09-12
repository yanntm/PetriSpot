#pragma once
#include <memory>
#include <numeric>
#include <optional>
#include <stdexcept>
#include "core/SparsePetriNet.h"

namespace petri::reduction {
/** Original guards and the free-distribution family represented by each
 * current coordinate. In-memory only; see algorithm.md. */
template<class T> class Enabling {
public:
  struct Origin {
    MatrixCol<T> pre;
    std::vector<std::string> names;
  };
  struct Group {
    size_t parent, size;
    std::optional<T> fixed;
  };
  std::shared_ptr<const Origin> origin;
  std::vector<Group> groups;
  std::vector<size_t> current;

  explicit Enabling(const SparsePetriNet<T>& net)
      : origin(std::make_shared<Origin>(Origin{net.getFlowPT(), net.getTnames()})),
        current(net.getPlaceCount()) {
    std::iota(current.begin(), current.end(), size_t{0});
    for (size_t p : current) groups.push_back({p, 1, std::nullopt});
  }
  size_t root(size_t p) const {
    while (groups.at(p).parent != p) p = groups[p].parent;
    return p;
  }
  void fuse(size_t p, size_t kept) {
    size_t a = root(current.at(p)), b = root(current.at(kept));
    if (a == b || groups[a].fixed || groups[b].fixed)
      throw std::logic_error("Invalid enabling-group fusion");
    if (groups[a].size > groups[b].size) std::swap(a, b);
    groups[a].parent = b;
    groups[b].size += groups[a].size;
    current[kept] = b;
  }
  void constant(size_t p, T marking) {
    groups[root(current.at(p))].fixed = marking;
  }
  Enabling compact(const std::vector<size_t>& map) const {
    Enabling out = *this;
    out.current.clear();
    for (size_t p = 0; p < map.size(); ++p) {
      if (map[p] == size_t(-1)) continue;
      if (out.current.size() <= map[p]) out.current.resize(map[p] + 1);
      out.current[map[p]] = root(current.at(p));
    }
    return out;
  }
};
}
