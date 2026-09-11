#pragma once
#include <algorithm>
#include "reduction/Workspace.h"
namespace petri::reduction {
/** Structural implicit place (ITS-Tools inducedBy sufficient condition).
 * p has one unit feeder h and one unit consumer f; h is a two-output fork and
 * f a two-input join. The other join input q starts empty and has a unit
 * single-feeder causal chain originating in h. Hence each token consumed at
 * q requires a prior h, and f cannot demand a p token that has not arrived.
 * p may initially have extra tokens, which only weaken its redundant guard.
 * Depth limits completeness, not the implication. No state/token counting. */
struct ImplicitForkJoin {
  static constexpr const char* name = "implicit-fork-join";
  template<class T> static bool induced(const Workspace<T>& w, size_t p, size_t cause,
                                       size_t depth, std::vector<size_t>& path) {
    if (!depth || w.marks[p] != 0 || std::find(path.begin(), path.end(), p) != path.end()) return false;
    const auto& in = w.producers.getColumn(p);
    if (in.size() != 1 || in.valueAt(0) != 1) return false;
    size_t feeder = in.keyAt(0);
    if (feeder == cause) return true;
    path.push_back(p);
    const auto& pre = w.pre.getColumn(feeder);
    bool found = false;
    for (size_t i = 0; i < pre.size(); ++i)
      if (pre.valueAt(i) == 1 && induced(w, pre.keyAt(i), cause, depth - 1, path)) { found = true; break; }
    path.pop_back();
    return found;
  }
  template<class T> static void apply(Workspace<T>& w) {
    if (w.config.goal == Goal::STATESPACE) return;
    std::vector<bool> removed(w.places.size());
    for (size_t p = w.places.size(); p-- > 0;) {
      if (p % 256 == 0 && w.stop()) return;
      if (!w.liveP[p] || w.observed[p]) continue;
      const auto& in = w.producers.getColumn(p); const auto& out = w.consumers.getColumn(p);
      if (in.size() != 1 || out.size() != 1 || in.valueAt(0) != 1 || out.valueAt(0) != 1) continue;
      size_t h = in.keyAt(0), f = out.keyAt(0);
      if (w.post.getColumn(h).size() != 2 || w.pre.getColumn(f).size() != 2) continue;
      const auto& join = w.pre.getColumn(f);
      size_t index = join.keyAt(0) == p ? 1 : 0;
      if (join.valueAt(index) != 1 || removed[join.keyAt(index)]) continue;
      std::vector<size_t> path;
      if (induced(w, join.keyAt(index), h, w.config.implicitDepth, path)) removed[p] = true;
    }
    for (size_t p = 0; p < removed.size(); ++p) if (removed[p]) w.retirePlace(p);
  }
};
}
