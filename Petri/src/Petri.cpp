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
#include "MatrixExporter.h"
#include "PNMLExport.h"  

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
const string MINFLOWS = "--minBasis";
const string EXPORT_MATRIX = "--exportAsMatrix";
const string NORMALIZE_PNML = "--normalizePNML";  

#define DEFAULT_TIMEOUT 150

#ifndef VAL
// default to 64 bit
#define VAL long
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
      << "  --minBasis           Force minimization of the semi flows basis (default: false).\n"
      << "  --exportAsMatrix=<file>  Export the incidence matrix to <file> in sparse format.\n"
      << "  --normalizePNML[=<file>] Exports a normalized PNML model to <file> (default: <input_path>.norm.pnml).\n"
      << "                       All places and transitions have id and name set to p0,p1... and t0,t1...\n"
      << "                       respectively, graphical and tool-specific information is removed.\n"
      << "  -q                   Quiet mode: Suppress detailed invariant output.\n"
      << "  -t <seconds>         Set timeout for computations (default: "
      << DEFAULT_TIMEOUT << "s).\n"
      << "  --noSingleSignRow    Disable single sign row heuristic in invariant computation.\n"
      << "  --noTrivialCull      Disable test for equal or empty columns before algorithm.\n"
      << "  --pivot=<strategy>   Set pivot strategy for elimination heuristic:\n"
      << "                       - best: Optimize for best pivot (default).\n"
      << "                       - worst: Use worst pivot (for testing).\n"
      << "                       - first: Use first valid pivot.\n"
      << "  --loopLimit=<n>      Limit elimination loops to <n> iterations (-1 for no limit).\n\n"
      << "Notes:\n" << "  - P-flows and P-semiflows are mutually exclusive.\n"
      << "  - T-flows and T-semiflows are mutually exclusive.\n"
      << "  - Invariant options enable invariant analysis.\n"
      << "  - Output files (e.g., .dot, .mat, .pnml) are written to the current directory.\n\n"
      << "Examples:\n" << "  petri -i model.pnml --draw\n"
      << "  petri -i model.pnml --findDeadlock -t 300\n"
      << "  petri -i model.pnml --Psemiflows -q --minBasis\n"
      << "  petri -i model.pnml --exportAsMatrix=model.mat\n"
      << "  petri -i x/y/model.pnml --normalizePNML\n"
      << "  petri -i model.pnml --normalizePNML=normalized.pnml\n";
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
  bool useCulling = true;
  bool minimizeFlows = false;
  std::string exportMatrixFile;
  std::string normalizePnmlFile; 
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
    } else if (std::string (argv[i]) == "--noTrivialCull") {
      useCulling = false;
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
    } else if (std::string (argv[i]) == MINFLOWS) {
      minimizeFlows = true;
    } else if (std::string (argv[i]).substr (0, 17) == "--exportAsMatrix=") {
      exportMatrixFile = std::string (argv[i]).substr (17);
    } else if (std::string (argv[i]).substr (0, 15) == "--normalizePNML") {
      std::string arg = std::string(argv[i]);
      if (arg == NORMALIZE_PNML) {
        // Default output file: append .norm.pnml to input path
        normalizePnmlFile = "";
      } else if (arg.substr(0, 16) == "--normalizePNML=") {
        normalizePnmlFile = arg.substr(16);
      }
    } else {
      std::cout << "[WARNING   ] Option : " << argv[i] << " not recognized"
          << std::endl;
    }
  }

  EliminationHeuristic heur (useSingleSignRow, pivotStrategy, loopLimit,
                             useCulling, minimizeFlows);

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

    // Handle PNML normalization and export
    if (!normalizePnmlFile.empty() || normalizePnmlFile == "") {
      // Set default output file if not specified
      std::string outputFile = normalizePnmlFile;
      if (outputFile.empty()) {
        outputFile = modelPath + ".norm.pnml";
      }
      pn->normalize();  // Call our normalize function
      PNMLExport<VAL>::transform(*pn, outputFile);
    }

    if (draw) {
      std::string title = "Petri Net: " + pn->getName ();
      std::string filename = FlowPrinter<VAL>::drawNet (
          *pn, title, std::set<size_t> (), std::set<size_t> (),
          std::numeric_limits<size_t>::max ());
      std::string targetFile = pn->getName () + ".dot";
      if (std::rename (filename.c_str (), targetFile.c_str ()) == 0) {
        std::cout << "Renamed output to " << targetFile << std::endl;
      } else {
        std::cerr << "Warning: Could not rename " << filename << " to "
            << targetFile << std::endl;
      }
    }

    // Export matrix if requested
    if (!exportMatrixFile.empty ()) {
      MatrixCol<VAL> sumMatrix = MatrixCol<VAL>::sumProd (-1, pn->getFlowPT (),
                                                          1, pn->getFlowTP ());
      MatrixExporter<VAL>::exportMatrix (sumMatrix, exportMatrixFile);
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
        auto invar = InvariantMiddle<VAL>::computePInvariants (sumMatrix,
                                                               psemiflows,
                                                               timeout, heur);
        std::cout << "Computed " << invar.getColumnCount () << " P "
            << (psemiflows ? "semi" : "") << "flows in "
            << std::chrono::duration_cast<std::chrono::milliseconds> (
                std::chrono::steady_clock::now () - time).count () << " ms."
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
        auto invar = InvariantMiddle<VAL>::computePInvariants (sumMatrix,
                                                               tsemiflows,
                                                               timeout, heur);
        std::cout << "Computed " << invar.getColumnCount () << " T "
            << (tsemiflows ? "semi" : "") << "flows in "
            << std::chrono::duration_cast<std::chrono::milliseconds> (
                std::chrono::steady_clock::now () - time).count () << " ms."
            << std::endl;
        if (!quiet) {
          std::vector<VAL> emptyVector;
          InvariantMiddle<VAL>::printInvariant (invar, pn->getTnames (),
                                                emptyVector);
        }
      }
    }

    delete pn;

  } catch (const char *e) {
    std::cout << e << std::endl;
    return 1;
  }

  std::cout << "Total runtime "
      << std::chrono::duration_cast<std::chrono::milliseconds> (
          std::chrono::steady_clock::now () - runtime).count () << " ms."
      << std::endl;

  return 0;
}
