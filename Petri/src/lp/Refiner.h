/*
 * Refiner.h
 *
 * The interface between a candidate solution and the refinements of the
 * state equation (traps, read arcs, predecessor, integrality: algorithm.md
 * section 4). A refiner examines a solved problem and either accepts the
 * solution, adds rows (a cut) or splits the problem in two. `refine` runs
 * the loop: solve, ask every refiner in turn, re-solve until acceptance,
 * infeasibility of every branch, or the budget. No refiner exists yet: the
 * loop with an empty list is a single solve.
 */
#ifndef PETRI_LP_REFINER_H_
#define PETRI_LP_REFINER_H_

#include <cstddef>
#include <memory>
#include <vector>

#include "lp/LpProblem.h"
#include "lp/Simplex.h"

namespace petri::lp
{

struct Verdict
{
  enum class Kind
  {
    Accept, Cut, Split
  };
  Kind kind = Kind::Accept;
  std::vector<Row> cuts;              // Cut: rows to add to the problem
  std::vector<std::vector<Row>> branches; // Split: one problem per branch, each with these rows added
};

class Refiner
{
public:
  virtual ~Refiner () = default;
  /** Judge a solution of p; the verdict says what to add before solving again. */
  virtual Verdict examine (const LpProblem &p, const LpResult &solution) = 0;
};

struct RefineOutcome
{
  LpResult result;      // the accepted solution, or the last status when none
  size_t solves = 0;    // linear programs solved
  size_t cuts = 0;
  size_t branches = 0;
  bool infeasible = false; // every branch infeasible: the goal is unreachable
};

/**
 * Depth-first over the branches, at most `maxSolves` linear programs. The
 * base problem is shared; a node of the tree is the list of rows added on
 * top of it (its overlay), so a branch costs its own rows and nothing of the
 * base.
 */
inline RefineOutcome refine (const LpProblem &base, const std::vector<std::unique_ptr<Refiner>> &refiners,
                             const LpLimits &limits, size_t maxSolves = 64, bool debug = false)
{
  RefineOutcome out;
  std::vector<std::vector<Row>> stack;
  stack.push_back ({});
  bool sawLimit = false;
  while (!stack.empty () && out.solves < maxSolves) {
    std::vector<Row> overlay = std::move (stack.back ());
    stack.pop_back ();
    Simplex simplex (limits);
    simplex.debug = debug;
    LpResult r = simplex.solve (base, overlay);
    ++out.solves;
    out.result = r;
    if (r.status == LpStatus::Infeasible) continue;
    if (!r.feasible ()) { sawLimit = true; continue; }
    bool accepted = true;
    for (const auto &ref : refiners) {
      Verdict v = ref->examine (base, r);
      if (v.kind == Verdict::Kind::Accept) continue;
      accepted = false;
      if (v.kind == Verdict::Kind::Cut) {
        for (auto &row : v.cuts) overlay.push_back (std::move (row));
        ++out.cuts;
        stack.push_back (std::move (overlay));
      } else {
        for (auto &rows : v.branches) {
          std::vector<Row> b = overlay;
          for (auto &row : rows) b.push_back (std::move (row));
          stack.push_back (std::move (b));
          ++out.branches;
        }
      }
      break;
    }
    if (accepted) return out;
  }
  out.infeasible = stack.empty () && !sawLimit && out.solves > 0;
  return out;
}

} // namespace petri::lp

#endif /* PETRI_LP_REFINER_H_ */
