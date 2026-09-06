/*
 * LpDriver.h
 *
 * --lp: every property of --props through the state equation (lp/). A goal
 * whose relaxation is infeasible in every conjunct is answered with a FORMULA
 * line (TECHNIQUES STATE_EQUATION); a feasible one yields the Parikh vector
 * of its shortest rational solution, written as (parikh NAME (t k)...) to
 * --lpHints for a later --hints. A bound property is maximised and its
 * relaxation optimum printed as an upper bound. One line per property
 * reports the status, size and cost of the solve; a summary closes.
 */
#ifndef PETRI_CLI_LPDRIVER_H_
#define PETRI_CLI_LPDRIVER_H_

#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "cli/Options.h"
#include "cli/WalkDriver.h"
#include "core/SparsePetriNet.h"
#include "expr/Property.h"
#include "lp/Parikh.h"
#include "lp/Refiner.h"
#include "lp/StateEquation.h"

namespace petri::cli
{

/** A transition name as an s-expression atom, quoted when it is not a plain symbol. */
inline std::string sexprAtom (const std::string &name)
{
  bool plain = !name.empty ();
  for (char c : name)
    if (!(std::isalnum (static_cast<unsigned char> (c)) || c == '_' || c == '-' || c == '.' || c == ':')) plain = false;
  if (plain) return name;
  std::string q = "\"";
  for (char c : name) {
    if (c == '"' || c == '\\') q += '\\';
    q += c;
  }
  return q + "\"";
}

template<typename T>
  void runLp (const Options &o, const SparsePetriNet<T> &pn)
  {
    using petri::expr::PropertyKind;
    using petri::lp::LpStatus;
    auto t0 = std::chrono::steady_clock::now ();
    std::vector<petri::expr::Property> props = loadProperties (o, pn);
    petri::lp::StateEquation<T> se (pn);
    std::cout << "State equation: " << pn.getPlaceCount () << " places, " << pn.getTransitionCount ()
        << " transitions, effects in " << millisSince (t0) << " ms." << std::endl;
    std::ofstream hints;
    if (!o.lpHintsFile.empty ()) hints.open (o.lpHintsFile);
    std::vector<std::unique_ptr<petri::lp::Refiner>> refiners; // none yet
    size_t proved = 0, hinted = 0, unknown = 0, skipped = 0, bounds = 0;
    size_t totalSolves = 0, totalPivots = 0;
    for (const auto &prop : props) {
      auto tp = std::chrono::steady_clock::now ();
      petri::lp::LpLimits limits;
      limits.hasDeadline = true;
      limits.deadline = tp + std::chrono::milliseconds (static_cast<long> (o.lpTime * 1000));
      if (prop.kind == PropertyKind::Unsupported || prop.kind == PropertyKind::Deadlock) {
        std::cout << "LP " << prop.name << ": skipped (" << to_string (prop.kind) << ")." << std::endl;
        ++skipped;
        continue;
      }
      if (prop.kind == PropertyKind::Bound) {
        petri::lp::LpProblem p = se.build ({});
        if (!se.objectiveMaximise (prop.boundForm (), p)) {
          std::cout << "LP " << prop.name << ": bound constant " << se.initialValue (prop.boundForm ()) << "." << std::endl;
          ++bounds;
          continue;
        }
        petri::lp::Simplex simplex (limits);
        petri::lp::LpResult r = simplex.solve (p);
        totalSolves++;
        totalPivots += r.pivots;
        std::cout << "LP " << prop.name << ": bound " << to_string (r.status);
        if (r.status == LpStatus::Optimal)
          std::cout << ", relaxation maximum " << (se.initialValue (prop.boundForm ()) - r.objective);
        std::cout << ", " << r.pivots << " pivots, " << millisSince (tp) << " ms." << std::endl;
        ++bounds;
        continue;
      }
      petri::expr::Expression goal = prop.goal ();
      std::vector<petri::lp::Conjunct> cases;
      if (!petri::lp::StateEquation<T>::toConjuncts (goal, cases)) {
        std::cout << "LP " << prop.name << ": " << (cases.size () > petri::lp::StateEquation<T>::MAX_CASES ? "declined, too many conjuncts" : "goal is constant false") << "." << std::endl;
        ++skipped;
        continue;
      }
      bool feasible = false, limited = false;
      size_t solves = 0, pivots = 0, constantFalse = 0;
      petri::expr::ParikhHint hint;
      for (const auto &cj : cases) {
        // an atom no transition moves is decided at the initial marking
        bool dead = false;
        for (const auto &a : cj)
          if (se.isConstant (a) && !petri::expr::compare (se.initialValue (a), a.op, a.constant)) dead = true;
        if (dead) { ++constantFalse; continue; }
        petri::lp::LpProblem p = se.build (cj);
        petri::lp::StateEquation<T>::objectiveShortest (p);
        petri::lp::RefineOutcome out = petri::lp::refine (std::move (p), refiners, limits);
        solves += out.solves;
        pivots += out.result.pivots;
        if (out.result.feasible ()) {
          feasible = true;
          hint = petri::lp::toParikh (out.result.x);
          break;
        }
        if (!out.infeasible) limited = true;
      }
      totalSolves += solves;
      totalPivots += pivots;
      std::cout << "LP " << prop.name << ": " << cases.size () << " conjunct" << (cases.size () > 1 ? "s" : "") << ", ";
      if (feasible) {
        ++hinted;
        std::cout << "feasible, Parikh of " << hint.counts.size () << " transitions, " << hint.total () << " firings";
        if (hints.is_open ()) {
          hints << "(parikh " << sexprAtom (prop.name);
          for (const auto &c : hint.counts) hints << " (" << sexprAtom (pn.getTnames ()[c.first]) << " " << c.second << ")";
          hints << ")\n";
        }
      } else if (limited) {
        ++unknown;
        std::cout << "unknown (budget)";
      } else {
        ++proved;
        std::cout << "infeasible";
        std::cout << std::endl << "FORMULA " << prop.name << " " << (prop.kind == PropertyKind::Invariant ? "TRUE" : "FALSE")
            << " TECHNIQUES STATE_EQUATION LP";
      }
      std::cout << ", " << solves << " solves, " << pivots << " pivots, " << millisSince (tp) << " ms." << std::endl;
    }
    std::cout << "LP summary: " << props.size () << " properties, " << proved << " proved unreachable, " << hinted
        << " with a Parikh vector, " << unknown << " unknown, " << bounds << " bounds, " << skipped << " skipped; "
        << totalSolves << " solves, " << totalPivots << " pivots, " << millisSince (t0) << " ms." << std::endl;
  }

} // namespace petri::cli

#endif /* PETRI_CLI_LPDRIVER_H_ */
