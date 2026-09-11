#pragma once
#include <array>
#include "reduction/rules/SourceTransition.h"
#include "reduction/rules/ConstantPlace.h"
#include "reduction/rules/DuplicateTransition.h"
#include "reduction/rules/DuplicatePlace.h"
#include "reduction/rules/EmptySiphon.h"
#include "reduction/rules/NoEffect.h"
#include "reduction/rules/SinkPlace.h"
#include "reduction/rules/PrefixOfInterest.h"
#include "reduction/rules/FreeSCC.h"
#include "reduction/rules/LoopBack.h"
#include "reduction/rules/TrivialPost.h"
#include "reduction/rules/PostAgglo.h"
#include "reduction/rules/PreAgglo.h"
#include "reduction/rules/ImplicitForkJoin.h"

#include "reduction/rules/FutureEquivalent.h"
#include "reduction/rules/FreeAgglo.h"
#include "reduction/rules/PartialFreeAgglo.h"
#include "reduction/rules/PartialPostAgglo.h"
#include "reduction/rules/RedundantComposition.h"
#include "reduction/rules/ScalarTransition.h"
#include "reduction/rules/SinkTransition.h"
#include "reduction/rules/InitialTokenMove.h"

namespace petri::reduction {
struct RuleStats { const char* name = ""; size_t passes = 0, edits = 0; };
struct NoTrace { static constexpr bool enabled = false; };

/** Java reachability/deadlock structural phase order. Free SCC and prefix lead;
 * local cleanup, prefix, implicit, SCC and trivial/simple post reach stability.
 * Pre/future/general/complex fallbacks precede siphons, redundancy, free/partial
 * agglomeration and token movement. Four growing outer rounds stop expansion.
 * STATESPACE runs its own short loop of record-maintaining rules instead.
 * Trace policies are compile-time optional: no event construction when off.
 * An enabled policy supplies before(rule,workspace) / after(rule,workspace).
 * The policy owns filtering, bounded capture and synchronous output. */
template<class T, class Trace = NoTrace> class Coordinator {
  Workspace<T>& w;
  Trace& trace;
  template<class Rule> void run(size_t index) {
    if (w.stop() || w.deadlock.has_value()) return;
    size_t before = w.changes;
    if constexpr (Trace::enabled) trace.before(Rule::name, w);
    Rule::apply(w);
    if constexpr (Trace::enabled) trace.after(Rule::name, w);
    auto& stat = stats[index]; stat.name = Rule::name; ++stat.passes;
    stat.edits += w.changes - before;
  }
public:
  std::array<RuleStats, 26> stats {};
  size_t passes = 0;
  Coordinator(Workspace<T>& net, Trace& observer) : w(net), trace(observer) {}
  void execute() {
    if (w.config.goal == Goal::NONE) return;
    if (w.config.goal == Goal::STATESPACE) {
      // the rules audited for a counting record (algorithm.md section 4), to stability
      size_t before;
      do {
        if (++passes > w.config.maxPasses || w.stop()) { w.limited = true; return; }
        before = w.changes;
        run<ConstantPlace>(0); run<DuplicateTransition>(1); run<NoEffect>(3); run<FreeSCC>(11);
      } while (w.changes != before);
      return;
    }
    run<SourceTransition>(25);
    size_t initial = w.changes;
    run<FreeSCC>(11);
    if (w.changes != initial) {
      run<NoEffect>(3); run<SinkTransition>(16); run<DuplicateTransition>(1); run<SourceTransition>(25); run<ScalarTransition>(17);
    }
    run<PrefixOfInterest>(5);
    size_t growing = 0;
    while (passes < w.config.maxPasses && !w.stop() && !w.deadlock.has_value()) {
      size_t countBefore = std::count(w.liveT.begin(), w.liveT.end(), true);
      size_t before;
      do {
        if (++passes > w.config.maxPasses || w.stop() || w.deadlock.has_value()) break;
        before = w.changes;
        run<ConstantPlace>(0); run<SinkPlace>(4); run<DuplicatePlace>(2); run<LoopBack>(12);
        run<NoEffect>(3); run<SinkTransition>(16); run<DuplicateTransition>(1); run<SourceTransition>(25); run<ScalarTransition>(17);
        run<PrefixOfInterest>(5); run<ImplicitForkJoin>(9);
        if (w.changes != before) run<FreeSCC>(11);
        size_t trivial = w.changes;
        run<TrivialPost>(7);
        if (w.changes == trivial) run<SimplePostAgglo>(13);
      } while (w.changes != before);
      if (passes > w.config.maxPasses || w.stop() || w.deadlock.has_value()) break;
      before = w.changes;
      run<PreAgglo>(10);
      bool rename = false;
      for (size_t t = 0; t < w.transitions.size(); ++t)
        if (w.liveT[t] && w.transitions[t].size() >= w.config.maxNameBytes) { rename = true; break; }
      if (rename) {
        size_t index = 0;
        for (size_t t = 0; t < w.transitions.size(); ++t)
          if (w.liveT[t]) w.transitions[t] = "t" + std::to_string(index++);
      }
      if (w.changes == before) run<FutureEquivalent>(18);
      if (w.changes == before) run<PostAgglo>(8);
      if (w.changes == before) run<ComplexPostAgglo>(14);
      if (w.changes == before) run<ComplexPreAgglo>(15);
      if (w.changes == before) run<FreeSCC>(11);
      if (w.changes == before) run<PrefixOfInterest>(5);
      run<EmptySiphon>(6); run<ConstantPlace>(0); run<SinkPlace>(4); run<DuplicatePlace>(2); run<LoopBack>(12);
      if (w.changes == before) run<RedundantComposition>(19);
      if (w.changes == before) run<FreeAgglo>(20);
      if (w.changes == before) run<ComplexFreeAgglo>(21);
      if (w.changes == before) run<PartialFreeAgglo>(22);
      if (w.changes == before && w.config.goal == Goal::REACHABILITY) run<PartialPostAgglo>(23);
      if (w.changes == before) {
        run<InitialTokenMove>(24); run<ConstantPlace>(0); run<SinkPlace>(4); run<DuplicatePlace>(2); run<LoopBack>(12);
      }
      if (w.changes == before) return;
      size_t countAfter = std::count(w.liveT.begin(), w.liveT.end(), true);
      growing = countAfter > countBefore ? growing + 1 : 0;
      if (growing > 3) { w.limited = true; return; }
    }
    if (passes >= w.config.maxPasses) w.limited = true;
  }
};
}
