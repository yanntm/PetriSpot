/*
 * Hunt.h
 *
 * The witness hunt for E[a U b] and E[a W b] from a state: guarded random
 * walks with restarts under a step budget, epsilon-greedy on the marking
 * distance to suf(b) when the formula has one; for W a deadlock or a lasso
 * closes the hunt. Sub-obligations (a temporal a or b) go back to the
 * checker. See algorithm.md.
 */
#ifndef PETRI_CTL_HUNT_H_
#define PETRI_CTL_HUNT_H_

#include <cstdint>
#include <random>
#include <unordered_map>
#include <vector>

#include "ctl/Cursor.h"
#include "ctl/Evidence.h"
#include "ctl/Verdict.h"
#include "expr/CtlFormula.h"
#include "expr/CtlSimplify.h"
#include "expr/Distance.h"

namespace petri::ctl
{

template<typename T>
  class Checker;

template<typename T>
  class Hunt
  {
    using CtlFormula = petri::expr::CtlFormula;
    using Expression = petri::expr::Expression;

    Checker<T> &chk;
    std::mt19937_64 &rng;

    /**
     * How many times t can fire in a row from the marking: it stays enabled as
     * long as every place it depletes keeps its input weight; 1 when it
     * depletes nothing.
     */
    static uint64_t maxFirings (const Cursor<T> &cur, uint32_t t)
    {
      const SparseArray<T> &eff = cur.net->effect (t);
      const SparseArray<T> &pre = cur.net->pre (t);
      uint64_t k = std::numeric_limits<uint64_t>::max ();
      for (size_t i = 0; i < eff.size (); ++i) {
        T d = eff.valueAt (i);
        if (d >= 0) continue;
        size_t p = eff.keyAt (i);
        T have = cur.marking.get (p) - pre.get (p);
        k = std::min (k, static_cast<uint64_t> (have / -d) + 1);
      }
      return k == std::numeric_limits<uint64_t>::max () ? 1 : k;
    }

    /** Epsilon-greedy choice: the sampled successor nearest to suf by the marking distance. */
    uint32_t choose (Cursor<T> &cur, const Expression *suf, const Budget &b)
    {
      size_t n = cur.enabled.size ();
      uint32_t pick = cur.enabled.at (rng () % n);
      if (!suf || rng () % 100 < b.epsilon) return pick;
      uint64_t best = std::numeric_limits<uint64_t>::max ();
      size_t tries = b.sample ? std::min (b.sample, n) : n;
      for (size_t i = 0; i < tries; ++i) {
        uint32_t t = b.sample && b.sample < n ? cur.enabled.at (rng () % n) : cur.enabled.at (i);
        uint64_t d = cur.marking.peek (cur.net->effect (t), [&] (const petri::walk::Marking<T> &m) {
          return petri::expr::distance (*suf, m);
        });
        if (d < best || (d == best && (rng () & 1))) {
          best = d;
          pick = t;
        }
      }
      return pick;
    }

  public:
    Hunt (Checker<T> &c, std::mt19937_64 &r)
        : chk (c), rng (r)
    {
    }

    /** A found path proves node at every state along it (the suffix is a witness): into the memo. */
    void rememberPath (const Cursor<T> &start, const std::vector<uint32_t> &path, const CtlFormula &node, Verdict found)
    {
      Cursor<T> re (start);
      for (uint32_t t : path) {
        chk.remember (re.state (), node, found);
        re.fire (t);
      }
    }

    /**
     * Hunt a witness of E[a U b] (weak: E[a W b]) from start. True with the
     * evidence, False when the start state decides it, Unknown when the
     * budget ends first. node is the formula the hunt decides and found its
     * verdict when a witness exists (False when hunting the negation of
     * node); every state on the witness path is remembered with it.
     */
    Verdict run (const Cursor<T> &start, const CtlFormula &a, const CtlFormula &b, bool weak, Evidence &ev,
                 const Budget &budget, const CtlFormula &node, Verdict found)
    {
      Stats &st = chk.stats ();
      ++st.hunts;
      Evidence sub;
      Verdict rb = chk.solve (start, b, &sub);
      if (rb == Verdict::True) {
        ev.what = "holds at the state itself";
        ev.kids.push_back (std::move (sub));
        return Verdict::True;
      }
      Verdict ra = chk.solve (start, a, nullptr);
      if (ra != Verdict::True) return ra == Verdict::False && rb == Verdict::False ? Verdict::False : Verdict::Unknown;
      if (start.deadlock ()) {
        if (weak) {
          ev.what = "EG: the state is a deadlock where the left side holds";
          return Verdict::True;
        }
        return rb == Verdict::False ? Verdict::False : Verdict::Unknown;
      }
      Expression sufExpr = petri::expr::ctlSuf (b);
      const Expression *suf = sufExpr.isConstant () ? nullptr : &sufExpr;
      // the run's state is local: a sub-obligation re-enters run
      std::vector<uint32_t> path;
      std::unordered_map<SparseArray<T>, size_t> seen; // W: states of the run, by their index on the path
      uint64_t steps = 0;
      uint64_t runs = 0;
      while (steps < budget.huntSteps && !chk.timedOut ()) {
        ++st.huntRuns;
        // saturation on every other run: the chosen transition is repeated while it stays enabled, every
        // intermediate state checked as any other (stacks of tokens move in a few choices)
        bool saturated = budget.saturate && (runs++ & 1);
        uint64_t repeat = 0;
        uint32_t last = 0;
        Cursor<T> cur (start);
        path.clear ();
        seen.clear ();
        if (weak) seen.emplace (cur.state (), 0);
        for (uint64_t len = 0; len < budget.runLength && steps < budget.huntSteps; ++len) {
          if (cur.deadlock ()) {
            if (weak) {
              ev.what = "EG: path to a deadlock, the left side holding all along";
              ev.path = path;
              rememberPath (start, path, node, found);
              return Verdict::True;
            }
            break;
          }
          uint32_t t;
          if (repeat > 0 && cur.enabled.isEnabled (last)) {
            t = last;
            --repeat;
          } else {
            t = choose (cur, suf, budget);
            if (saturated) {
              repeat = maxFirings (cur, t) - 1;
              last = t;
            }
          }
          cur.fire (t);
          path.push_back (t);
          ++steps;
          ++st.huntSteps;
          Evidence kid;
          if (chk.solve (cur, b, &kid) == Verdict::True) {
            ev.what = weak ? "EW: path to the right side" : "EU: path to the right side";
            ev.path = path;
            ev.kids.push_back (std::move (kid));
            rememberPath (start, path, node, found);
            return Verdict::True;
          }
          if (chk.solve (cur, a, nullptr) != Verdict::True) break;
          if (weak) {
            auto ins = seen.emplace (cur.state (), path.size ());
            if (!ins.second) {
              ev.what = "EG: lasso, the left side holding all along";
              ev.path = path;
              ev.loopIndex = static_cast<long> (ins.first->second);
              rememberPath (start, path, node, found);
              return Verdict::True;
            }
          }
        }
      }
      return Verdict::Unknown;
    }
  };

} // namespace petri::ctl

#endif /* PETRI_CTL_HUNT_H_ */
