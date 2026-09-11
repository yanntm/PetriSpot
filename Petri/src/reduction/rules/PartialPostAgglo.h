#pragma once
#include "reduction/Composition.h"
namespace petri::reduction {
/** Java rulePartialPostAgglo. Empty unobserved p controls all consumers by a
 * single unit guard; both visible and invisible consumers exist. With disjoint
 * unit feeders, compose just invisible consumers. Feeders, visible consumers
 * and p remain. The reachability schedule restricts p to at most one feeder;
 * Java's direct DEADLOCK invocation allows several, but reduce() does not call
 * this rule for that goal. SI scheduling is deferred. */
struct PartialPostAgglo {
  static constexpr const char* name = "partial-post-agglomeration";
  template<class T> static void apply(Workspace<T>& w) {
    if ((w.config.goal != Goal::REACHABILITY && w.config.goal != Goal::DEADLOCK) || !w.config.agglomeration) return;
    std::vector<size_t> candidates;
    for (size_t p = 0; p < w.places.size(); ++p) if (w.liveP[p] && !w.observed[p] && w.marks[p] == 0) {
      const auto& row = w.consumers.getColumn(p);
      bool ok = true, visible = false, invisible = false;
      for (size_t i = 0; i < row.size(); ++i) {
        size_t t = row.keyAt(i);
        if (w.pre.getColumn(t).size() != 1 || row.valueAt(i) != 1) { ok = false; break; }
        if (w.visible(t)) visible = true; else invisible = true;
      }
      if (ok && visible && invisible) candidates.push_back(p);
    }
    std::vector<size_t> removed;
    size_t added = 0;
    for (size_t p : candidates) {
      if (w.stop()) break;
      auto in = w.producers.getColumn(p), out = w.consumers.getColumn(p);
      if (w.config.goal != Goal::DEADLOCK && in.size() > 1) continue;
      bool ok = true;
      for (size_t i = 0; i < in.size(); ++i)
        if (in.valueAt(i) != 1 || out.get(in.keyAt(i))) { ok = false; break; }
      if (!ok) continue;
      for (size_t i = 0; i < out.size(); ++i) {
        size_t f = out.keyAt(i);
        if (w.visible(f) || w.pre.getColumn(f).size() != 1 || out.valueAt(i) != 1) continue;
        for (size_t j = 0; j < in.size(); ++j) {
          size_t h = in.keyAt(j);
          auto pre = composeColumn(w.pre.getColumn(h), w.pre.getColumn(f), T{1}, p);
          auto post = composeColumn(w.post.getColumn(h), w.post.getColumn(f), T{1}, p);
          auto name = w.composedName(h, f);
          w.appendTransition(std::move(pre), std::move(post), std::move(name)); ++added;
        }
        removed.push_back(f);
      }
    }
    if (added) for (size_t t : removed) w.retireTransition(t);
  }
};
}
