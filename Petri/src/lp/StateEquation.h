/*
 * StateEquation.h
 *
 * From a net and a goal to linear programs over the transition counts. The
 * marking is the affine form s_p(x) = m0_p + sum_t C[p,t] x_t, never a
 * variable: an atom over places becomes one row over transitions, the
 * non-negativity of every touched place one row each. A goal in disjunctive
 * normal form gives one problem per conjunct; strict comparisons are
 * tightened to the integers before they enter (s < k is s <= k - 1). See
 * algorithm.md sections 1 and 2.
 */
#ifndef PETRI_LP_STATEEQUATION_H_
#define PETRI_LP_STATEEQUATION_H_

#include <cstddef>
#include <cstdint>
#include <map>
#include <vector>

#include "core/MatrixCol.h"
#include "core/SparseArray.h"
#include "core/SparsePetriNet.h"
#include "expr/Expression.h"
#include "lp/LpProblem.h"

namespace petri::lp
{

/** A conjunct of the goal: atoms with the comparison already oriented (negations pushed in). */
using Conjunct = std::vector<petri::expr::LinearAtom>;

template<typename T>
  class StateEquation
  {
    const SparsePetriNet<T> &net;
    MatrixCol<T> effects;   // by transition: post - pre
    MatrixCol<T> effectsT;  // by place: the effects of every transition on it
    std::vector<T> m0;

  public:
    /** Conjuncts beyond this many and the goal is declined rather than enumerated. */
    static constexpr size_t MAX_CASES = 64;

    explicit StateEquation (const SparsePetriNet<T> &n)
        : net (n), m0 (n.getMarks ())
    {
      size_t nt = n.getTransitionCount ();
      effects = MatrixCol<T> (n.getPlaceCount (), 0);
      effects.reserveColumns (nt);
      for (size_t t = 0; t < nt; ++t)
        effects.appendColumn (SparseArray<T>::sumProd (1, n.getFlowTP ().getColumn (t), -1, n.getFlowPT ().getColumn (t)));
      effectsT = effects.transpose ();
    }

    const MatrixCol<T>& effectMatrix () const
    {
      return effects;
    }
    size_t transitionCount () const
    {
      return net.getTransitionCount ();
    }

    /**
     * The goal as a disjunction of conjuncts of oriented atoms; false when the
     * goal is constant false or has more than MAX_CASES conjuncts. A constant
     * true goal yields one empty conjunct.
     */
    static bool toConjuncts (const petri::expr::Expression &goal, std::vector<Conjunct> &out)
    {
      out.clear ();
      return dnf (goal, false, out) && out.size () <= MAX_CASES;
    }

    /** The linear program of one conjunct: non-negativity of the marking plus the atoms; no objective. */
    LpProblem build (const Conjunct &atoms) const
    {
      LpProblem p (net.getTransitionCount ());
      for (size_t pl = 0; pl < effectsT.getColumnCount (); ++pl) {
        const SparseArray<T> &line = effectsT.getColumn (pl);
        if (line.size () == 0) continue;
        Row r;
        r.coeffs.reserve (line.size ());
        for (size_t i = 0; i < line.size (); ++i)
          r.coeffs.emplace_back (static_cast<uint32_t> (line.keyAt (i)), static_cast<long long> (line.valueAt (i)));
        r.lo = -static_cast<double> (m0[pl]);
        r.hi = INF;
        p.addRow (std::move (r));
      }
      for (const auto &a : atoms) {
        Row r;
        long long rhs;
        if (!atomRow (a, r, rhs)) continue; // an atom over untouched places only: constant, checked by the caller
        setRange (a.op, rhs, r);
        p.addRow (std::move (r));
      }
      return p;
    }

    /** Minimise the total firing count: a short Parikh vector. */
    static void objectiveShortest (LpProblem &p)
    {
      p.setObjective (std::vector<double> (p.columns, 1.0));
    }

    /** Maximise a linear form over places, as a form over transitions; false when no transition moves it. */
    bool objectiveMaximise (const petri::expr::LinearAtom &form, LpProblem &p) const
    {
      Row r;
      long long shift;
      if (!atomRow (form, r, shift)) return false;
      std::vector<double> c (p.columns, 0.0);
      for (const auto &e : r.coeffs) c[e.first] = static_cast<double> (e.second);
      p.setMaximise (c);
      return true;
    }

    /** Value of a form at the initial marking, the constant part of its affine expression. */
    long long initialValue (const petri::expr::LinearAtom &form) const
    {
      long long v = 0;
      for (const auto &t : form.terms) v += t.second * static_cast<long long> (m0[t.first]);
      return v;
    }

    /** Whether an atom holds or fails regardless of the transitions (no transition touches its places). */
    bool isConstant (const petri::expr::LinearAtom &a) const
    {
      for (const auto &t : a.terms)
        if (effectsT.getColumn (t.first).size () > 0) return false;
      return true;
    }

  private:
    /**
     * The row of an atom over transitions: coefficient of x_t is sum_p a_p C[p,t];
     * rhs is the constant minus the form's value at m0. False when no
     * transition touches the atom's places.
     */
    bool atomRow (const petri::expr::LinearAtom &a, Row &r, long long &rhs) const
    {
      std::map<uint32_t, long long> acc;
      for (const auto &term : a.terms) {
        const SparseArray<T> &line = effectsT.getColumn (term.first);
        for (size_t i = 0; i < line.size (); ++i)
          acc[static_cast<uint32_t> (line.keyAt (i))] += term.second * static_cast<long long> (line.valueAt (i));
      }
      r.coeffs.clear ();
      for (const auto &e : acc)
        if (e.second != 0) r.coeffs.push_back (e);
      rhs = a.constant - initialValue (a);
      return !r.coeffs.empty ();
    }

    /** lo/hi of a row from the comparison, tightened to the integers. NEQ never reaches here (split by dnf). */
    static void setRange (petri::expr::Cmp op, long long rhs, Row &r)
    {
      using petri::expr::Cmp;
      double v = static_cast<double> (rhs);
      switch (op) {
      case Cmp::GE: r.lo = v; r.hi = INF; break;
      case Cmp::GT: r.lo = v + 1; r.hi = INF; break;
      case Cmp::LE: r.lo = -INF; r.hi = v; break;
      case Cmp::LT: r.lo = -INF; r.hi = v - 1; break;
      case Cmp::EQ: r.lo = v; r.hi = v; break;
      case Cmp::NEQ: r.lo = -INF; r.hi = INF; break;
      }
    }

    /** Disjunctive normal form with negations pushed to the atoms; a NEQ atom splits in two. */
    static bool dnf (const petri::expr::Expression &e, bool neg, std::vector<Conjunct> &out)
    {
      using petri::expr::Expression;
      using petri::expr::Cmp;
      switch (e.kind) {
      case Expression::Kind::True:
      case Expression::Kind::False: {
        bool val = (e.kind == Expression::Kind::True) != neg;
        if (val) out.push_back ({});
        return val;
      }
      case Expression::Kind::Not:
        return dnf (e.children[0], !neg, out);
      case Expression::Kind::Atom: {
        petri::expr::LinearAtom a = e.atom;
        if (neg) a.op = petri::expr::negate (a.op);
        if (a.op == Cmp::NEQ) {
          petri::expr::LinearAtom lt = a, gt = a;
          lt.op = Cmp::LT;
          gt.op = Cmp::GT;
          out.push_back ({ lt });
          out.push_back ({ gt });
        } else {
          out.push_back ({ a });
        }
        return true;
      }
      case Expression::Kind::Or:
      case Expression::Kind::And: {
        bool isOr = (e.kind == Expression::Kind::Or) != neg; // De Morgan under negation
        if (isOr) {
          bool any = false;
          for (const auto &c : e.children) {
            std::vector<Conjunct> sub;
            if (dnf (c, neg, sub)) {
              any = true;
              for (auto &cj : sub) out.push_back (std::move (cj));
              if (out.size () > MAX_CASES) return true; // the caller sees the overflow
            }
          }
          return any;
        }
        // a conjunction: the product of the children's conjuncts
        std::vector<Conjunct> acc { {} };
        for (const auto &c : e.children) {
          std::vector<Conjunct> sub;
          if (!dnf (c, neg, sub)) return false;
          std::vector<Conjunct> next;
          for (const auto &left : acc)
            for (const auto &right : sub) {
              Conjunct cj = left;
              cj.insert (cj.end (), right.begin (), right.end ());
              next.push_back (std::move (cj));
              if (next.size () > MAX_CASES) { out = std::move (next); return true; }
            }
          acc = std::move (next);
        }
        for (auto &cj : acc) out.push_back (std::move (cj));
        return true;
      }
      }
      return false;
    }
  };

} // namespace petri::lp

#endif /* PETRI_LP_STATEEQUATION_H_ */
