#pragma once
#include "reduction/Composition.h"
namespace petri::reduction {
/** Java rulePostAgglo F-continuation. Every F is controlled solely by p;
 * H/F are disjoint and weights permit whole continuations without splitting
 * tokens between choices. At most one side is visible; reachability additionally
 * restricts effects of H when F is visible, temporal modes require invisible F.
 * Marked p is allowed in the complex phase with one invisible continuation.
 * Simple/non-complex/complex are scheduling choices, not different theorems. */
struct PostAgglo {
  static constexpr const char* name = "post-agglomeration";
  template<class T> static void apply(Workspace<T>& w, bool complex = false, bool simple = false) {
    if (!permitsAgglomeration(w.config)) return;
    size_t applications = 0;
    for (size_t p = 0; p < w.places.size(); ++p) {
      if (w.stop()) return;
      if (!w.liveP[p] || w.observed[p]) continue;
      const auto& in = w.producers.getColumn(p); const auto& out = w.consumers.getColumn(p);
      if ((in.size() == 0) || (out.size() == 0)) continue;
      if (simple && (in.size() != 1 || out.size() != 1 || w.post.getColumn(in.keyAt(0)).size() > 1)) continue;
      if (!complex && in.size() > 1 && out.size() > 1) continue;
      if (w.marks[p] != 0 && (out.size() > 1 || !complex)) continue;
      if (w.config.goal == Goal::SI_CTL && out.size() > 1) continue;
      if (in.size() > 1 && out.size() > 1
          && (in.size() > (w.config.postCrossProductLimit - 1) / out.size())) continue;
      bool ok = true, visibleH = false, visibleF = false;
      std::vector<size_t> hs, fs;
      for (size_t i = 0; i < out.size(); ++i) {
        size_t f = out.keyAt(i);
        if (w.pre.getColumn(f).size() != 1 || w.marks[p] % out.valueAt(i) != 0) { ok = false; break; }
        visibleF |= w.visible(f); fs.push_back(f);
      }
      for (size_t i = 0; ok && i < in.size(); ++i) {
        size_t h = in.keyAt(i);
        if (w.pre.getColumn(h).get(p)) { ok = false; break; }
        for (size_t j = 0; j < out.size(); ++j) {
          T fed = in.valueAt(i), consumed = out.valueAt(j);
          if (fed % consumed != 0 || (fed > consumed && out.size() > 1)) { ok = false; break; }
        }
        visibleH |= w.visible(h); hs.push_back(h);
      }
      if (!ok || (visibleH && visibleF)) continue;
      if (w.marks[p] != 0 && visibleF) continue;
      Goal goal = w.config.goal;
      if ((goal == Goal::SI_LTL || goal == Goal::LI_LTL || goal == Goal::SI_CTL) && visibleF) continue;
      if (goal == Goal::REACHABILITY && visibleF) {
        for (size_t h : hs) {
          size_t effects = 0;
          const auto& pre = w.pre.getColumn(h); const auto& post = w.post.getColumn(h);
          for (size_t i = 0; i < pre.size(); ++i) if (pre.valueAt(i) != post.get(pre.keyAt(i))) ++effects;
          for (size_t i = 0; i < post.size(); ++i) if (pre.get(post.keyAt(i)) == 0) ++effects;
          if (effects > 1) { ok = false; break; }
        }
      }
      if (ok && agglomerate(w, p, hs, fs) && complex && ++applications >= w.config.complexPostApplications) return;
    }
  }
};
struct SimplePostAgglo {
  static constexpr const char* name = "post-agglomeration-simple";
  template<class T> static void apply(Workspace<T>& w) { PostAgglo::apply(w, false, true); }
};
struct ComplexPostAgglo {
  static constexpr const char* name = "post-agglomeration-complex";
  template<class T> static void apply(Workspace<T>& w) { PostAgglo::apply(w, true); }
};
}
