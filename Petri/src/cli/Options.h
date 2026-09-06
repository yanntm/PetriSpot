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
  std::string netFile;       // --net, PNET binary net
  std::string loadKERSFile;  // --loadKERS, incidence matrix instead of a net
  bool quiet = false;
  int timeout = 150;         // seconds, invariants and single walks

  // exports and inspection
  bool draw = false;
  std::string exportMatrixFile;
  std::string exportAsKERSFile;
  std::string exportNetFile;  // --exportNet, PNET
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
  std::string propsSyntax = "auto"; // auto (by extension) | mcc | sexpr
  std::string hintsFile;            // --hints, s-expression (parikh NAME (t k)...) forms
  std::string hintStrategies = "parikh,parikh,relaxed,bestfirst"; // pool for a focus with a hint
  long query = -1;
  bool printProps = false;
  std::string printPropsFormat = "infix"; // infix | sexpr | sexpr-index
  bool findDeadlock = false;
  petri::walk::WalkBudget budget;
  uint64_t seed = static_cast<uint64_t> (std::chrono::steady_clock::now ().time_since_epoch ().count ());
  bool trace = false;
  bool printUnknown = false;
  std::string strategy = "random";
  std::string strategies;
  unsigned threads = 1;
  unsigned epsilon = 10;
  size_t sample = 0;
  uint64_t stall = 0;
  uint64_t debugSteps = 0;
  size_t share = 32;
  unsigned shareProb = 50;
  long totalTime = 0;
  bool escalate = false;   // raise the step budget rather than concede a round
  long roundTime = 1;
  long sweepTime = 1; // multi-target round before the focused rounds (0: none)
  std::string sweepChoice = "rare"; // transition choice of the sweep: rare (rarity and age) or random
  uint64_t runTime = 1000;   // wall clock budget of one run before a restart, ms (0: none)
  uint64_t noveltyStall = 100000; // sweep: restart after this many steps without a new transition (0: never)
  size_t partition = 0; // sweeps over at least this many targets are split between threads (0: never)

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
  auto *pnml = in->add_option ("-i", o.modelPath, "Input model file (PNML, ISO/IEC 15909-2).");
  in->add_option ("--net", o.netFile, "Input net in PNET binary format (places p<i>, transitions t<i>; see INTEROP.md).")
      ->excludes (pnml);
  in->add_option ("--loadKERS", o.loadKERSFile,
                  "Sparse integer matrix in KERS format instead of a net, for --Pflows/--Psemiflows (rows are the "
                  "variables) or --Tflows/--Tsemiflows (transposed internally).");
  in->add_flag ("-q", o.quiet, "Quiet: no invariant listing, fewer walk reports.");
  in->add_option ("-t", o.timeout, "Timeout in seconds for a computation (default 150).");

  auto *ex = app.add_option_group ("Exports and inspection");
  ex->add_flag ("--draw", o.draw, "Write the net as <model>.dot (Graphviz).");
  ex->add_option ("--exportAsMatrix", o.exportMatrixFile, "Write the incidence matrix in ASCII sparse format.");
  ex->add_option ("--exportAsKERS", o.exportAsKERSFile, "Write the incidence matrix in KERS format.");
  ex->add_option ("--exportNet", o.exportNetFile, "Write the net in PNET binary format.");
  ex->add_option_function<std::vector<std::string>> ("--normalizePNML",
      [&o] (const std::vector<std::string> &v) {
        o.normalizePNML = true;
        if (!v.empty () && !v[0].empty ()) o.normalizePNMLFile = v[0];
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
  wk->add_option ("--props", o.propsFile,
                  "Property file: MCC XML (.xml) or s-expressions (any other extension, see INTEROP.md).");
  wk->add_option ("--propsSyntax", o.propsSyntax, "Property syntax: auto (by extension, default), mcc, sexpr.")
      ->check (CLI::IsMember ({ "auto", "mcc", "sexpr" }));
  wk->add_option ("--hints", o.hintsFile,
                  "Hints file: (parikh NAME (t k)...) forms naming properties of --props (see INTEROP.md).");
  wk->add_option ("--hintStrategies", o.hintStrategies,
                  "Pool for a focus that has a hint (default parikh,parikh,relaxed,bestfirst).");
  wk->add_option ("--query", o.query, "Select the n-th property of the file (0-based); default all.");
  wk->add_option_function<std::vector<std::string>> ("--printProps",
      [&o] (const std::vector<std::string> &v) {
        o.printProps = true;
        if (!v.empty () && !v[0].empty ()) o.printPropsFormat = v[0];
      }, "Print the selected properties and exit: infix (parsed and normalised, default), sexpr (names), "
         "sexpr-index (p<i>/t<i>).")
      ->expected (0, 1)->check (CLI::IsMember ({ "infix", "sexpr", "sexpr-index" }));
  wk->add_flag ("--findDeadlock", o.findDeadlock, "Look for a deadlock.");
  wk->add_option ("--walkSteps", o.budget.maxSteps, "Step budget of a walk (default: until timeout).");
  wk->add_option ("--runLength", o.budget.runLength, "Steps between restarts from the initial marking (default 1000000).");
  wk->add_option ("--seed", o.seed, "Random seed (default: clock).");
  wk->add_flag ("--trace", o.trace, "Record the witness trace, verify it by replay and print it as a WITNESS line.");
  wk->add_flag ("--printUnknown", o.printUnknown, "At exit, print an UNKNOWN line for every property left without verdict.");
  wk->add_option ("--strategy", o.strategy,
                  "random (default) | rare | bestfirst | structural | relaxed | parikh | deadlock | sync, each optionally with a +sat suffix "
                  "(fire the chosen transition as many times as the marking allows).")
      ->check (CLI::Validator ([] (std::string &s) {
        std::string base = s;
        if (base.size () > 4 && base.compare (base.size () - 4, 4, "+sat") == 0) base.erase (base.size () - 4);
        for (const char *k : { "random", "rare", "bestfirst", "structural", "relaxed", "parikh", "deadlock", "sync" })
          if (base == k) return std::string ();
        return "unknown strategy " + s;
      }, "STRATEGY"));
  wk->add_option ("--threads", o.threads, "Parallel walkers (default 1); first witness wins.");
  wk->add_option ("--strategies", o.strategies,
                  "Pool assigned round-robin to threads: name[:epsilon[:stall]],... (default with several "
                  "threads: random,bestfirst,structural,relaxed).");
  wk->add_option ("--share", o.share, "Shared pool of up to n promising markings for restarts (default 32, 0 disables).");
  wk->add_option ("--shareProb", o.shareProb, "Percentage of restarts drawn from the pool (default 50).");
  wk->add_option ("--totalTime", o.totalTime,
                  "Global budget in seconds for all properties, walked in rounds of growing per-property budget.");
  wk->add_option ("--roundTime", o.roundTime, "Per-property budget of the first round, seconds (default 1).");
  wk->add_flag ("--escalate", o.escalate,
                "When a round solves nothing and every walk ended on its step budget, multiply the step budget "
                "and walk on, instead of concluding that more time cannot help.");
  wk->add_option ("--sweepTime", o.sweepTime,
                  "With at least two open properties, a first round of random walks checking all of them at once, "
                  "seconds (default 1, 0 disables).");
  wk->add_option ("--sweepChoice", o.sweepChoice,
                  "Transition choice of the sweep: rare (the least fired of a few sampled, the oldest enabled "
                  "among equals) or random (default rare).");
  wk->add_option ("--runTime", o.runTime, "Wall clock budget of one run before a restart, ms (default 1000, 0 none).");
  wk->add_option ("--noveltyStall", o.noveltyStall,
                  "Sweep: restart after this many steps without firing a new transition (default 100000, 0 never).");
  wk->add_option ("--partition", o.partition,
                  "Sweeps over at least this many properties are split between the threads, each checking its "
                  "share: fewer checks per step, but a thread walks past the finds another owns (default 0, never).");
  wk->add_option ("--epsilon", o.epsilon, "Heuristic strategies: percentage of random moves (default 10).");
  wk->add_option ("--sample", o.sample, "bestfirst: candidates scored per step (default all).");
  wk->add_option ("--stall", o.stall, "Heuristic strategies: restart after n steps without improvement.");
  wk->add_option ("--debugSteps", o.debugSteps, "Trace the first n relaxed-plan decisions on stderr.");
}

} // namespace petri::cli

#endif /* PETRI_CLI_OPTIONS_H_ */
