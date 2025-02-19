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

#define DEFAULT_TIMEOUT 150

#ifndef VAL
// default to 32 bit
#define VAL int
#endif

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
  int timeout = DEFAULT_TIMEOUT;
  // New elimination heuristic configuration parameters.
  bool useSingleSignRow = true;
  EliminationHeuristic::PivotStrategy pivotStrategy = EliminationHeuristic::PivotStrategy::FindBest;
  ssize_t loopLimit = -1;

  if (argc == 1 || argc > 6) {
    cerr << "usage: petri -i model.pnml [options]\n";
    exit (1);
  }

  for (int i = 1; i < argc; i++) {
    if (argv[i] == PATH) {
      modelPath = argv[++i];
    } else if (argv[i] == TIMEOUT) {
       timeout = atoi(argv[++i]);
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
    } else if (std::string(argv[i]) == "--noSingleSignRow") {
      useSingleSignRow = false;
    } else if (std::string(argv[i]).substr(0, 8) == "--pivot=") {
      std::string pivotStr = std::string(argv[i]).substr(8);
      if (pivotStr == "best") {
        pivotStrategy = EliminationHeuristic::PivotStrategy::FindBest;
      } else if (pivotStr == "worst") {
        pivotStrategy = EliminationHeuristic::PivotStrategy::FindWorst;
      } else if (pivotStr == "first") {
        pivotStrategy = EliminationHeuristic::PivotStrategy::FindFirst;
      } else {
        std::cerr << "Unknown pivot strategy: " << pivotStr << std::endl;
        exit(1);
      }
    } else if (std::string(argv[i]).substr(0, 12) == "--loopLimit=") {
      std::string limitStr = std::string(argv[i]).substr(12);
      try {
        loopLimit = std::stoll(limitStr);
      } catch (const std::exception &e) {
        std::cerr << "Invalid loopLimit value: " << limitStr << std::endl;
        exit(1);
      }
    } else {
      std::cout << "[WARNING   ] Option : " << argv[i] << " not recognized"
          << std::endl;
    }
  }

  // Create the elimination heuristic configuration.
    EliminationHeuristic heur(useSingleSignRow, pivotStrategy, loopLimit);

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

  if(modelPath.empty()){
      std::cerr << "Error: no model file specified." << std::endl;
      return 1;
  }

  std::ifstream file(modelPath);
  if (!file.good()) {
      std::cerr << "Error: file not found: " << modelPath << std::endl;
      return 1;
  }


  try {
    SparsePetriNet<VAL> *pn = loadXML<VAL> (modelPath);

// 		std::cout << "PN : " ;
// 		std::cout << "\nPre : " << std::endl;
// 		pn->getFlowPT().print(std::cout);
// 		std::cout << "\nPost : " << std::endl;
// 		pn->getFlowTP().print(std::cout);
// 		std::cout << std::endl;

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
        MatrixCol<VAL> sumMatrix =  MatrixCol<VAL>::sumProd(-1, pn->getFlowPT(), 1, pn->getFlowTP());
        unordered_set<SparseArray<VAL>> invar = InvariantMiddle<VAL>::computePInvariants (sumMatrix, psemiflows , timeout, heur);

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

        MatrixCol<VAL> sumMatrix =  MatrixCol<VAL>::sumProd(-1, pn->getFlowPT(), 1, pn->getFlowTP()).transpose();

        unordered_set<SparseArray<VAL>> invarT = InvariantMiddle<VAL>::computePInvariants (sumMatrix, tsemiflows, timeout);

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
