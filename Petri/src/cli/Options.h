/*
 * Options.h
 *
 * Every command line setting of the petri binaries, with its default, and the
 * CLI11 declaration of the options. See README.md.
 */
#ifndef PETRI_CLI_OPTIONS_H_
#define PETRI_CLI_OPTIONS_H_

#include <chrono>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "cli11/CLI11.hpp"
#include "invariants/Heuristic.h"
#include "walk/Walker.h"

namespace petri::cli
{

struct Options
{
  // input
  std::string modelPath;     // -i, PNML
  std::string loadKERSFile;  // --loadKERS, incidence matrix instead of a net
  bool quiet = false;
  int timeout = 150;         // seconds, invariants and single walks

  // exports and inspection
  bool draw = false;
  std::string exportMatrixFile;
  std::string exportAsKERSFile;
  bool normalizePNML = false;
  std::string normalizePNMLFile; // empty: <model>.norm.pnml
  bool netStats = false;

  // invariants
  bool pflows = false;
  bool psemiflows = false;
  bool tflows = false;
  bool tsemiflows = false;
  bool useSingleSignRow = true;
  bool useCulling = true;
  bool minimizeFlows = false;
  bool useQPlusBasis = false;
  bool useCompression = false;
  EliminationHeuristic::PivotStrategy pivot = EliminationHeuristic::PivotStrategy::FindBest;
  ssize_t loopLimit = 500;
  std::string basisKERSFile;

  // reachability walk
  std::string propsFile;
  long query = -1;
  bool printProps = false;
  bool findDeadlock = false;
  petri::walk::WalkBudget budget;
  uint64_t seed = static_cast<uint64_t> (std::chrono::steady_clock::now ().time_since_epoch ().count ());
  bool trace = false;
  std::string strategy = "random";
  std::string strategies;
  unsigned threads = 1;
  unsigned epsilon = 10;
  size_t sample = 0;
  uint64_t stall = 0;
  uint64_t debugSteps = 0;
  size_t share = 0;
  unsigned shareProb = 50;
  long totalTime = 0;
  long roundTime = 1;

  bool invariants () const
  {
    return pflows || psemiflows || tflows || tsemiflows;
  }

  EliminationHeuristic heuristic (bool compression) const
  {
    return EliminationHeuristic (useSingleSignRow, pivot, loopLimit, useCulling,
                                 minimizeFlows, useQPlusBasis, compression);
  }
};

/** Declare every option of o on app. Values land in o when app.parse runs. */
inline void addOptions (CLI::App &app, Options &o)
{
  app.set_help_flag ("-h,--help", "Print this help and exit.");
  app.footer ("Examples:\n"
              "  petri -i model.pnml --Psemiflows -q --minBasis\n"
              "  petri -i model.pnml --exportAsKERS=model.kers\n"
              "  petri --loadKERS=m.kers --Pflows --basisKERS=basis.kers\n"
              "  petri -i model.pnml --props=ReachabilityFireability.xml --threads=4 --totalTime=60\n"
              "  petri -i model.pnml --findDeadlock -t 30");

  auto *in = app.add_option_group ("Input");
  in->add_option ("-i", o.modelPath, "Input model file (PNML, ISO/IEC 15909-2).");
  in->add_option ("--loadKERS", o.loadKERSFile,
                  "Sparse integer matrix in KERS format instead of a net, for --Pflows/--Psemiflows (rows are the "
                  "variables) or --Tflows/--Tsemiflows (transposed internally).");
  in->add_flag ("-q", o.quiet, "Quiet: no invariant listing, fewer walk reports.");
  in->add_option ("-t", o.timeout, "Timeout in seconds for a computation (default 150).");

  auto *ex = app.add_option_group ("Exports and inspection");
  ex->add_flag ("--draw", o.draw, "Write the net as <model>.dot (Graphviz).");
  ex->add_option ("--exportAsMatrix", o.exportMatrixFile, "Write the incidence matrix in ASCII sparse format.");
  ex->add_option ("--exportAsKERS", o.exportAsKERSFile, "Write the incidence matrix in KERS format.");
  ex->add_option_function<std::vector<std::string>> ("--normalizePNML",
      [&o] (const std::vector<std::string> &v) {
        o.normalizePNML = true;
        if (!v.empty ()) o.normalizePNMLFile = v[0];
      }, "Write a normalised PNML (ids p0,p1... t0,t1..., no graphics); default <model>.norm.pnml.")
      ->expected (0, 1);
  ex->add_flag ("--netStats", o.netStats, "Print structural histograms of the net (arities, fan-out).");

  auto *inv = app.add_option_group ("Invariants");
  inv->add_flag ("--Pflows", o.pflows, "Generative basis of P-flows.");
  inv->add_flag ("--Psemiflows", o.psemiflows, "Generative basis of P-semiflows.");
  inv->add_flag ("--Tflows", o.tflows, "Generative basis of T-flows.");
  inv->add_flag ("--Tsemiflows", o.tsemiflows, "Generative basis of T-semiflows.");
  inv->add_flag ("--minBasis", o.minimizeFlows, "Minimise the semi-flow basis.");
  inv->add_option ("--basisKERS", o.basisKERSFile,
                   "Write the basis in KERS format (program-to-program mode: no listing).");
  inv->add_flag ("--useQPlusBasis", o.useQPlusBasis, "Q+ basis for semi-flows.");
  inv->add_flag ("--useCompression", o.useCompression, "Compress the semi-flow basis by permutations.");
  inv->add_flag ("!--noSingleSignRow", o.useSingleSignRow, "Disable the single sign row heuristic.");
  inv->add_flag ("!--noTrivialCull", o.useCulling, "Disable the removal of equal or empty columns.");
  std::map<std::string, EliminationHeuristic::PivotStrategy> pivots {
      { "best", EliminationHeuristic::PivotStrategy::FindBest },
      { "worst", EliminationHeuristic::PivotStrategy::FindWorst },
      { "first", EliminationHeuristic::PivotStrategy::FindFirst } };
  inv->add_option ("--pivot", o.pivot, "Pivot strategy of the elimination: best (default), worst, first.")
      ->transform (CLI::CheckedTransformer (pivots, CLI::ignore_case));
  inv->add_option ("--loopLimit", o.loopLimit, "Elimination loop limit (default 500, -1 unlimited).");

  auto *wk = app.add_option_group ("Reachability walk");
  wk->add_option ("--props", o.propsFile, "Property file: MCC XML (Reachability*.xml).");
  wk->add_option ("--query", o.query, "Select the n-th property of the file (0-based); default all.");
  wk->add_flag ("--printProps", o.printProps, "Print the selected properties (parsed, normalised) and exit.");
  wk->add_flag ("--findDeadlock", o.findDeadlock, "Look for a deadlock.");
  wk->add_option ("--walkSteps", o.budget.maxSteps, "Step budget of a walk (default: until timeout).");
  wk->add_option ("--runLength", o.budget.runLength, "Steps between restarts from the initial marking (default 1000000).");
  wk->add_option ("--seed", o.seed, "Random seed (default: clock).");
  wk->add_flag ("--trace", o.trace, "Record the witness trace, verify it by replay and print it.");
  wk->add_option ("--strategy", o.strategy, "random (default) | bestfirst | structural | relaxed.")
      ->check (CLI::IsMember ({ "random", "bestfirst", "structural", "relaxed" }));
  wk->add_option ("--threads", o.threads, "Parallel walkers (default 1); first witness wins.");
  wk->add_option ("--strategies", o.strategies,
                  "Pool assigned round-robin to threads: name[:epsilon[:stall]],... (default with several "
                  "threads: random,bestfirst,structural,relaxed).");
  wk->add_option ("--share", o.share, "Shared pool of up to n promising markings for restarts (default 0).");
  wk->add_option ("--shareProb", o.shareProb, "Percentage of restarts drawn from the pool (default 50).");
  wk->add_option ("--totalTime", o.totalTime,
                  "Global budget in seconds for all properties, walked in rounds of growing per-property budget.");
  wk->add_option ("--roundTime", o.roundTime, "Per-property budget of the first round, seconds (default 1).");
  wk->add_option ("--epsilon", o.epsilon, "Heuristic strategies: percentage of random moves (default 10).");
  wk->add_option ("--sample", o.sample, "bestfirst: candidates scored per step (default all).");
  wk->add_option ("--stall", o.stall, "Heuristic strategies: restart after n steps without improvement.");
  wk->add_option ("--debugSteps", o.debugSteps, "Trace the first n relaxed-plan decisions on stderr.");
}

} // namespace petri::cli

#endif /* PETRI_CLI_OPTIONS_H_ */
