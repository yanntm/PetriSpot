/*
 * Checker.h
 *
 * Local CTL checking over configurations (marking, formula node): a memo of
 * verdicts, EX / AX by successor enumeration, E nodes by hunts and A nodes
 * by regions, each falling back to the dual method on the node's negation.
 * Single-threaded. See algorithm.md.
 */
#ifndef PETRI_CTL_CHECKER_H_
#define PETRI_CTL_CHECKER_H_

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <memory>
#include <random>
#include <unordered_map>
#include <vector>

#include "ctl/Cursor.h"
#include "ctl/Evidence.h"
#include "ctl/Hunt.h"
#include "ctl/Region.h"
#include "ctl/Verdict.h"
#include "expr/CtlFormula.h"
#include "expr/CtlSimplify.h"
#include "walk/WalkNet.h"

namespace petri::ctl
{

template<typename T>
  class Checker
  {
    using CtlFormula = petri::expr::CtlFormula;
    using CtlOp = petri::expr::CtlOp;

    struct Key
    {
      SparseArray<T> marking;
      const CtlFormula *node;
      bool operator== (const Key &o) const
      {
        return node == o.node && marking == o.marking;
      }
    };
    struct KeyHash
    {
      size_t operator() (const Key &k) const
      {
        return k.marking.hash () * 31 + std::hash<const CtlFormula*> () (k.node);
      }
    };
    struct Entry
    {
      Verdict verdict = Verdict::Unknown;
      unsigned level = 0; // the round at which an unknown was last tried
    };

    const petri::walk::WalkNet<T> &net;
    std::mt19937_64 rng;
    Budget budget;
    Stats st;
    std::unordered_map<Key, Entry, KeyHash> memo;
    std::unordered_map<const CtlFormula*, std::unique_ptr<CtlFormula>> duals;
    Hunt<T> hunt;
    Region<T> region;
    unsigned depth = 0; // temporal nodes being solved above the current one
    const CtlFormula TRUE_F = CtlFormula::constant (true);
    const CtlFormula FALSE_F = CtlFormula::constant (false);

    /** The negation of a normal-form node, built once. */
    const CtlFormula& dualOf (const CtlFormula &f)
    {
      auto it = duals.find (&f);
      if (it == duals.end ()) it = duals.emplace (&f, std::make_unique<CtlFormula> (petri::expr::ctlDual (f))).first;
      return *it->second;
    }

    /** The until form of a temporal node: EF b is E[true U b], EG a is E[a W false]. */
    void untilForm (const CtlFormula &f, const CtlFormula *&a, const CtlFormula *&b, bool &weak)
    {
      switch (f.op) {
      case CtlOp::EF: case CtlOp::AF: a = &TRUE_F; b = &f.kids[0]; weak = false; break;
      case CtlOp::EG: case CtlOp::AG: a = &f.kids[0]; b = &FALSE_F; weak = true; break;
      case CtlOp::EU: case CtlOp::AU: a = &f.kids[0]; b = &f.kids[1]; weak = false; break;
      default: a = &f.kids[0]; b = &f.kids[1]; weak = true; break;
      }
    }

    Verdict solveNext (const Cursor<T> &s, const CtlFormula &f, Evidence *ev)
    {
      const CtlFormula &g = f.kids[0];
      bool exist = f.op == CtlOp::EX;
      if (s.deadlock ()) {
        // no successor: EX fails, AX holds vacuously (the MCC reading, as its-ctl and TAPAAL compute it)
        if (ev) ev->what = exist ? "EX: the state is a deadlock, no successor" : "AX: the state is a deadlock, no successor";
        return exist ? Verdict::False : Verdict::True;
      }
      std::vector<uint32_t> succ = s.enabled.transitions ();
      std::shuffle (succ.begin (), succ.end (), rng);
      bool unknown = false;
      for (uint32_t t : succ) {
        Cursor<T> c (s);
        c.fire (t);
        Evidence kid;
        Verdict v = solve (c, g, &kid);
        if (v == Verdict::Unknown) {
          unknown = true;
        } else if ((v == Verdict::True) == exist) {
          if (ev) {
            ev->what = exist ? "EX: a successor satisfies the operand" : "AX: a successor violates the operand";
            ev->path = { t };
            ev->kids.push_back (std::move (kid));
          }
          return exist ? Verdict::True : Verdict::False;
        }
      }
      if (unknown) return Verdict::Unknown;
      if (ev) ev->what = exist ? "EX: every successor violates the operand" : "AX: every successor satisfies the operand";
      return exist ? Verdict::False : Verdict::True;
    }

    Verdict solveBool (const Cursor<T> &s, const CtlFormula &f, Evidence *ev)
    {
      bool isAnd = f.op == CtlOp::And;
      // state children first: they are cheap and may decide the node
      std::vector<const CtlFormula*> order;
      for (const auto &k : f.kids) if (k.isState ()) order.push_back (&k);
      for (const auto &k : f.kids) if (!k.isState ()) order.push_back (&k);
      bool unknown = false;
      std::vector<Evidence> kids;
      for (const CtlFormula *k : order) {
        Evidence kid;
        Verdict v = solve (s, *k, &kid);
        if (v == Verdict::Unknown) {
          unknown = true;
        } else if ((v == Verdict::True) != isAnd) {
          // a false child of an and, a true child of an or: decided
          if (ev) {
            ev->what = isAnd ? "and: a child fails" : "or: a child holds";
            ev->kids.push_back (std::move (kid));
          }
          return isAnd ? Verdict::False : Verdict::True;
        } else if (!k->isState ()) {
          kids.push_back (std::move (kid));
        }
      }
      if (unknown) return Verdict::Unknown;
      if (ev) {
        ev->what = isAnd ? "and: every child holds" : "or: every child fails";
        ev->kids = std::move (kids);
      }
      return isAnd ? Verdict::True : Verdict::False;
    }

    /**
     * The budget of a hunt or a region at the current depth: the root budget
     * at the root, the nested budget for an obligation opened from a state
     * met along a hunt or in a region (the probe of one state).
     */
    Budget effectiveBudget () const
    {
      if (depth == 0) return budget;
      Budget b = budget;
      b.huntSteps = budget.nestedSteps;
      b.regionStates = budget.nestedStates;
      return b;
    }

    /** Depth accounting around a hunt or a region: what they open is nested. */
    struct Nested
    {
      unsigned &d;
      explicit Nested (unsigned &depth) : d (depth) { ++d; }
      ~Nested () { --d; }
    };

    Verdict runHunt (const Cursor<T> &s, const CtlFormula &a, const CtlFormula &b, bool weak, Evidence &ev,
                     const CtlFormula &node, Verdict found)
    {
      Budget eff = effectiveBudget ();
      Nested n (depth);
      return hunt.run (s, a, b, weak, ev, eff, node, found);
    }

    Verdict runRegion (const Cursor<T> &s, const CtlFormula &a, const CtlFormula &b, bool until, Evidence &ev,
                       const CtlFormula &node, Verdict proved)
    {
      Budget eff = effectiveBudget ();
      Nested n (depth);
      return region.run (s, a, b, until, ev, eff, node, proved);
    }

    /** The hunt for the negation of an A node f: a witness of the dual is a counter-example of f. */
    Verdict huntDual (const Cursor<T> &s, const CtlFormula &f, Evidence *ev)
    {
      const CtlFormula &d = dualOf (f);
      if (!petri::expr::isTemporal (d.op)) return Verdict::Unknown;
      const CtlFormula *da, *db;
      bool dweak;
      untilForm (d, da, db, dweak);
      Evidence viaDual;
      Verdict dv = runHunt (s, *da, *db, dweak, viaDual, f, Verdict::False);
      if (dv != Verdict::Unknown && ev) {
        ev->what = "by the negation, " + viaDual.what;
        ev->path = std::move (viaDual.path);
        ev->loopIndex = viaDual.loopIndex;
        ev->kids = std::move (viaDual.kids);
      }
      return negate (dv);
    }

    /** The region for the negation of an E node f: the dual proved is f refuted. */
    Verdict regionDual (const Cursor<T> &s, const CtlFormula &f, Evidence *ev)
    {
      const CtlFormula &d = dualOf (f);
      if (!petri::expr::isTemporal (d.op)) return Verdict::Unknown;
      const CtlFormula *da, *db;
      bool dweak;
      untilForm (d, da, db, dweak);
      Evidence viaDual;
      Verdict dv = runRegion (s, *da, *db, !dweak, viaDual, f, Verdict::False);
      if (dv != Verdict::Unknown && ev) {
        ev->what = "by the negation, " + viaDual.what;
        ev->path = std::move (viaDual.path);
        ev->loopIndex = viaDual.loopIndex;
        ev->kids = std::move (viaDual.kids);
      }
      return negate (dv);
    }

    /**
     * An E node: the hunt, then the region of its negation. An A node: the
     * cheap refutation first, a hunt for its negation, then the region.
     * Without allowDual (solving the negation of a node) only the direct
     * method runs, so that the two methods run once each.
     */
    Verdict solveTemporal (const Cursor<T> &s, const CtlFormula &f, Evidence *ev, bool allowDual)
    {
      const CtlFormula *a, *b;
      bool weak;
      untilForm (f, a, b, weak);
      bool exist = petri::expr::isExistential (f.op);
      Evidence direct;
      Verdict v = Verdict::Unknown;
      if (exist) {
        v = runHunt (s, *a, *b, weak, direct, f, Verdict::True);
        if (v != Verdict::Unknown) {
          if (ev) *ev = std::move (direct);
          return v;
        }
        if (!allowDual || timedOut ()) return Verdict::Unknown;
        return regionDual (s, f, ev);
      }
      if (allowDual) {
        v = huntDual (s, f, ev);
        if (v != Verdict::Unknown || timedOut ()) return v;
      }
      v = runRegion (s, *a, *b, !weak, direct, f, Verdict::True);
      if (v != Verdict::Unknown && ev) *ev = std::move (direct);
      return v;
    }

  public:
    Checker (const petri::walk::WalkNet<T> &n, uint64_t seed)
        : net (n), rng (seed), hunt (*this, rng), region (*this)
    {
    }

    void setBudget (const Budget &b)
    {
      budget = b;
    }
    const Budget& getBudget () const
    {
      return budget;
    }
    Stats& stats ()
    {
      return st;
    }
    bool timedOut () const
    {
      return std::chrono::steady_clock::now () >= budget.deadline;
    }
    size_t memoSize () const
    {
      return memo.size ();
    }

    /** A verdict known for (marking, node) without a solve: a state on a witness path or in a closed region. */
    void remember (const SparseArray<T> &m, const CtlFormula &node, Verdict v)
    {
      Entry &e = memo[Key { m, &node }];
      if (e.verdict == Verdict::Unknown) {
        e.verdict = v;
        ++st.propagated;
      }
    }

    /**
     * The verdict of f (normal form) at s. Evidence, when asked, describes
     * the witness or the counter-example. allowDual is false when solving the
     * negation of a node, so that the two methods run once each.
     */
    Verdict solve (const Cursor<T> &s, const CtlFormula &f, Evidence *ev, bool allowDual = true)
    {
      if (f.isState ()) {
        bool holds = f.evalState (s.marking, s.deadlock ());
        if (ev) ev->what = holds ? "state predicate holds" : "state predicate fails";
        return holds ? Verdict::True : Verdict::False;
      }
      if (f.op == CtlOp::Not) return negate (solve (s, f.kids[0], ev, allowDual));
      if (timedOut ()) return Verdict::Unknown;
      ++st.solves;
      Key key { s.state (), &f };
      auto it = memo.find (key);
      if (it != memo.end ()) {
        if (it->second.verdict != Verdict::Unknown) {
          ++st.memoHits;
          if (ev) ev->what = std::string ("(decided earlier: ") + to_string (it->second.verdict) + ")";
          return it->second.verdict;
        }
        if (it->second.level >= budget.level) {
          ++st.memoHits;
          return Verdict::Unknown;
        }
      } else {
        ++st.configurations;
      }
      Verdict v;
      if (f.op == CtlOp::And || f.op == CtlOp::Or) v = solveBool (s, f, ev);
      else if (f.op == CtlOp::EX || f.op == CtlOp::AX) v = solveNext (s, f, ev);
      else v = solveTemporal (s, f, ev, allowDual);
      Entry &e = memo[key];
      e.verdict = v;
      e.level = budget.level;
      return v;
    }
  };

} // namespace petri::ctl

#endif /* PETRI_CTL_CHECKER_H_ */
