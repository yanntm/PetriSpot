#include <iostream>
#include <string>
#include "Walker.h"
#include "PTNetLoader.h"
#include "InvariantMiddle.h"
#include <vector>
#include <fstream>
#include <unordered_set>
#include <chrono>
#include "Heuristic.h"
#include "FlowPrinter.h"

using namespace std;
using namespace petri;

const string FINDDEADLOCK = "--findDeadlock";
const string PFLOW = "--Pflows";
const string PSEMIFLOW = "--Psemiflows";
const string TFLOW = "--Tflows";
const string TSEMIFLOW = "--Tsemiflows";
const string PATH = "-i";
const string QUIET = "-q";
const string TIMEOUT = "-t";
const string DRAW = "--draw";

#define DEFAULT_TIMEOUT 150

#ifndef VAL
// default to 32 bit
#define VAL int
#endif

// Usage function to display help message
void usage ()
{
  std::cerr << "Usage: petri -i <model.pnml> [options]\n\n"
      << "PetriSpot: A tool for analyzing Petri nets from PNML files.\n\n"
      << "Required Arguments:\n"
      << "  -i <path>            Specify the input PNML model file.\n\n"
      << "Options:\n"
      << "  --draw               Export the Petri net to a <modelName>.dot file with no size limit.\n"
      << "  --findDeadlock       Run deadlock detection with forward and backward walks (up to 1M steps).\n"
      << "  --Pflows             Compute minimal P-flows (place invariants) of the net.\n"
      << "  --Psemiflows         Compute minimal P-semiflows (positive place invariants).\n"
      << "  --Tflows             Compute minimal T-flows (transition invariants).\n"
      << "  --Tsemiflows         Compute minimal T-semiflows (positive transition invariants).\n"
      << "  -q                   Quiet mode: Suppress detailed invariant output.\n"
      << "  -t <seconds>         Set timeout for computations (default: "
      << DEFAULT_TIMEOUT << "s).\n"
      << "  --noSingleSignRow    Disable single sign row heuristic in invariant computation.\n"
      << "  --pivot=<strategy>   Set pivot strategy for elimination heuristic:\n"
      << "                       - best: Optimize for best pivot (default).\n"
      << "                       - worst: Use worst pivot (for testing).\n"
      << "                       - first: Use first valid pivot.\n"
      << "  --loopLimit=<n>      Limit elimination loops to <n> iterations (-1 for no limit).\n\n"
      << "Notes:\n" << "  - P-flows and P-semiflows are mutually exclusive.\n"
      << "  - T-flows and T-semiflows are mutually exclusive.\n"
      << "  - Invariant options (--Pflows, --Psemiflows, --Tflows, --Tsemiflows) enable invariant analysis.\n"
      << "  - Output files (e.g., .dot) are written to the current working directory.\n\n"
      << "Examples:\n" << "  petri -i model.pnml --draw\n"
      << "  petri -i model.pnml --findDeadlock -t 300\n"
      << "  petri -i model.pnml --Psemiflows -q --pivot=first\n";
}

int main (int argc, char *argv[])
{
  std::string logMessage = "Running PetriSpot with arguments : [";
  int i = 0;
  for (i = 1; i < argc - 1; ++i) {
    logMessage += std::string (argv[i]) + ", ";
  }
  logMessage += (argc == 1 ? "Empty" : std::string (argv[i])) + "]";
  InvariantMiddle<VAL>::writeToLog (logMessage);
  auto runtime = std::chrono::steady_clock::now ();
  std::string modelPath;
  bool findDeadlock = false;
  bool pflows = false;
  bool tflows = false;
  bool psemiflows = false;
  bool tsemiflows = false;
  bool invariants = false;
  bool quiet = false;
  bool draw = false;
  int timeout = DEFAULT_TIMEOUT;
  bool useSingleSignRow = true;
  EliminationHeuristic::PivotStrategy pivotStrategy =
      EliminationHeuristic::PivotStrategy::FindBest;
  ssize_t loopLimit = -1;

  if (argc == 1) {
    usage ();
    exit (1);
  }

  for (int i = 1; i < argc; i++) {
    if (argv[i] == PATH) {
      modelPath = argv[++i];
    } else if (argv[i] == TIMEOUT) {
      timeout = atoi (argv[++i]);
    } else if (argv[i] == QUIET) {
      quiet = true;
    } else if (argv[i] == FINDDEADLOCK) {
      findDeadlock = true;
    } else if (argv[i] == PFLOW) {
      pflows = true;
      invariants = true;
    } else if (argv[i] == PSEMIFLOW) {
      psemiflows = true;
      invariants = true;
    } else if (argv[i] == TFLOW) {
      tflows = true;
      invariants = true;
    } else if (argv[i] == TSEMIFLOW) {
      tsemiflows = true;
      invariants = true;
    } else if (argv[i] == DRAW) {
      draw = true;
    } else if (std::string (argv[i]) == "--noSingleSignRow") {
      useSingleSignRow = false;
    } else if (std::string (argv[i]).substr (0, 8) == "--pivot=") {
      std::string pivotStr = std::string (argv[i]).substr (8);
      if (pivotStr == "best") {
        pivotStrategy = EliminationHeuristic::PivotStrategy::FindBest;
      } else if (pivotStr == "worst") {
        pivotStrategy = EliminationHeuristic::PivotStrategy::FindWorst;
      } else if (pivotStr == "first") {
        pivotStrategy = EliminationHeuristic::PivotStrategy::FindFirst;
      } else {
        std::cerr << "Unknown pivot strategy: " << pivotStr << std::endl;
        exit (1);
      }
    } else if (std::string (argv[i]).substr (0, 12) == "--loopLimit=") {
      std::string limitStr = std::string (argv[i]).substr (12);
      try {
        loopLimit = std::stoll (limitStr);
      } catch (const std::exception &e) {
        std::cerr << "Invalid loopLimit value: " << limitStr << std::endl;
        exit (1);
      }
    } else {
      std::cout << "[WARNING   ] Option : " << argv[i] << " not recognized"
          << std::endl;
    }
  }

  EliminationHeuristic heur (useSingleSignRow, pivotStrategy, loopLimit);

  if (pflows && psemiflows) {
    std::cout << "Cannot compute P flows and P semi-flows at the same time."
        << std::endl;
    return 1;
  }
  if (tflows && tsemiflows) {
    std::cout << "Cannot compute T flows and T semi-flows at the same time."
        << std::endl;
    return 1;
  }

  if (modelPath.empty ()) {
    std::cerr << "Error: no model file specified." << std::endl;
    usage ();
    return 1;
  }

  std::ifstream file (modelPath);
  if (!file.good ()) {
    std::cerr << "Error: file not found: " << modelPath << std::endl;
    return 1;
  }

  try {
    SparsePetriNet<VAL> *pn = loadXML<VAL> (modelPath);

    if (draw) {
      std::string title = "Petri Net: " + pn->getName ();
      // Updated call: std::set<int> -> std::set<size_t>, INT_MAX -> std::numeric_limits<size_t>::max()
      std::string filename = FlowPrinter<VAL>::drawNet (*pn,
                                                        title,
                                                        std::set<size_t> (), // Empty set for hlPlaces
          std::set<size_t> (), // Empty set for hlTrans
          std::numeric_limits < size_t > ::max () // Max size_t for no practical limit
              );
      std::string targetFile = pn->getName () + ".dot";
      if (std::rename (filename.c_str (), targetFile.c_str ()) == 0) {
        std::cout << "Renamed output to " << targetFile << std::endl;
      } else {
        std::cerr << "Warning: Could not rename " << filename << " to "
            << targetFile << std::endl;
      }
    }

    if (findDeadlock) {
      Walker<VAL> walk (*pn);
      if (walk.runDeadlockDetection (1000000, true, timeout)) {
        std::cout << "Deadlock found !" << std::endl;
        delete pn;
        return 0;
      } else {
        std::cout << "No deadlock found !" << std::endl;
      }
      if (walk.runDeadlockDetection (1000000, false, timeout)) {
        std::cout << "Deadlock found !" << std::endl;
      } else {
        std::cout << "No deadlock found !" << std::endl;
      }
    }

    if (invariants) {
      vector<int> repr;
      if (pflows || psemiflows) {
        auto time = std::chrono::steady_clock::now ();
        MatrixCol<VAL> sumMatrix = MatrixCol<VAL>::sumProd (-1,
                                                            pn->getFlowPT (), 1,
                                                            pn->getFlowTP ());
        unordered_set<SparseArray<VAL>> invar =
            InvariantMiddle<VAL>::computePInvariants (sumMatrix, psemiflows,
                                                      timeout, heur);

        std::cout << "Computed " << invar.size () << " P "
            << (psemiflows ? "semi" : "") << "flows in "
            << std::chrono::duration_cast < std::chrono::milliseconds
            > (std::chrono::steady_clock::now () - time).count () << " ms."
                << std::endl;
        if (!quiet) {
          InvariantMiddle<VAL>::printInvariant (invar, pn->getPnames (),
                                                (*pn).getMarks ());
        }
      }
      if (tflows || tsemiflows) {
        auto time = std::chrono::steady_clock::now ();
        MatrixCol<VAL> sumMatrix =
            MatrixCol<VAL>::sumProd (-1, pn->getFlowPT (), 1, pn->getFlowTP ()).transpose ();
        unordered_set<SparseArray<VAL>> invarT =
            InvariantMiddle<VAL>::computePInvariants (sumMatrix, tsemiflows,
                                                      timeout);

        std::cout << "Computed " << invarT.size () << " T "
            << (tsemiflows ? "semi" : "") << "flows in "
            << std::chrono::duration_cast < std::chrono::milliseconds
            > (std::chrono::steady_clock::now () - time).count () << " ms."
                << std::endl;
        if (!quiet) {
          std::vector<VAL> emptyVector;
          InvariantMiddle<VAL>::printInvariant (invarT, pn->getTnames (),
                                                emptyVector);
        }
      }
    }

    delete pn;

  } catch (const char *e) {
    std::cout << e << std::endl;
    return 1;
  }

  std::cout << "Total runtime " << std::chrono::duration_cast
      < std::chrono::milliseconds
      > (std::chrono::steady_clock::now () - runtime).count () << " ms."
          << std::endl;

  return 0;
}
