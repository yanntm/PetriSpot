#pragma once
#include <array>
#include "reduction/rules/ConstantPlace.h"
#include "reduction/rules/DuplicateTransition.h"
#include "reduction/rules/DuplicatePlace.h"
#include "reduction/rules/EmptySiphon.h"
#include "reduction/rules/NoEffect.h"
#include "reduction/rules/SinkPlace.h"
#include "reduction/rules/ReachabilityRelevance.h"
#include "reduction/rules/TrivialPost.h"
#include "reduction/rules/PostAgglo.h"
#include "reduction/rules/PreAgglo.h"
#include "reduction/rules/ImplicitForkJoin.h"

namespace petri::reduction {
struct RuleStats { const char* name = ""; size_t passes = 0, edits = 0; };
struct NoTrace { static constexpr bool enabled = false; };

/** Coordinator for the implemented subset. Cheap cleanup reaches stability
 * before siphons and agglomeration; each successful fallback returns to cleanup.
 * Trace policies are compile-time optional: no event construction when off.
 * An enabled policy supplies before(rule,workspace) / after(rule,workspace).
 * The policy owns filtering, bounded capture and synchronous output. */
template<class T, class Trace = NoTrace> class Coordinator {
  Workspace<T>& w;
  Trace& trace;
  template<class Rule> void run(size_t index) {
    if (w.stop()) return;
    size_t before = w.changes;
    if constexpr (Trace::enabled) trace.before(Rule::name, w);
    Rule::apply(w);
    if constexpr (Trace::enabled) trace.after(Rule::name, w);
    auto& stat = stats[index]; stat.name = Rule::name; ++stat.passes;
    stat.edits += w.changes - before;
  }
public:
  std::array<RuleStats, 11> stats {};
  size_t passes = 0;
  Coordinator(Workspace<T>& net, Trace& observer) : w(net), trace(observer) {}
  void execute() {
    if (w.config.goal == Goal::NONE) return;
    for (; passes < w.config.maxPasses && !w.stop(); ++passes) {
      size_t before = w.changes;
      run<ConstantPlace>(0); run<DuplicateTransition>(1); run<DuplicatePlace>(2);
      run<NoEffect>(3); run<SinkPlace>(4); run<ReachabilityRelevance>(5);
      if (w.changes != before) continue;
      run<EmptySiphon>(6);
      run<ImplicitForkJoin>(9);
      if (w.changes != before) continue;
      run<TrivialPost>(7);
      if (w.changes != before) continue;
      run<PostAgglo>(8);
      if (w.changes != before) continue;
      run<PreAgglo>(10);
      if (w.changes == before) return;
    }
    if (passes == w.config.maxPasses) w.limited = true;
  }
};
}
