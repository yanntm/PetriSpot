#include <iostream>
#include <string>
#include "walk/Walker.h"
#include "walk/Strategy.h"
#include "walk/Target.h"
#include "parse/PTNetLoader.h"
#include "invariants/InvariantMiddle.h"
#include <vector>
#include <fstream>
#include <unordered_set>
#include <chrono>
#include "invariants/Heuristic.h"
#include "io/FlowPrinter.h"
#include "io/MatrixExporter.h"
#include "io/PNMLExport.h"
#include "io/SparseMatrixIO.h"
#include "parse/mcc/PropertyLoader.h"

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
const string USEQPLUS = "--useQPlusBasis";
const string USECOMPRESSION = "--useCompression";
const string LOAD_KERS = "--loadKERS";
const string EXPORT_AS_KERS = "--exportAsKERS";
const string BASIS_KERS = "--basisKERS";
const string PROPS = "--props=";
const string QUERY = "--query=";
const string PRINT_PROPS = "--printProps";
const string WALK_STEPS = "--walkSteps=";
const string RUN_LENGTH = "--runLength=";
const string SEED = "--seed=";
const string TRACE = "--trace";

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
      << "  --exportAsMatrix=<file>  Export the incidence matrix to <file> in ASCII sparse format.\n"
      << "  --exportAsKERS=<file>    Export the incidence matrix to <file> in KERS binary format.\n"
      << "  --loadKERS=<file>        Load a sparse integer matrix in KERS binary format instead of PNML.\n"
      << "                           Use with --Pflows/--Psemiflows (rows=places, cols=transitions)\n"
      << "                           or --Tflows/--Tsemiflows (matrix is transposed internally).\n"
      << "                           Implies program-to-program mode when combined with --basisKERS.\n"
      << "  --basisKERS=<file>       Export computed invariant basis to <file> in KERS binary format.\n"
      << "                           Suppresses human-readable invariant output (program-to-program mode).\n"
      << "  --normalizePNML[=<file>] Exports a normalized PNML model to <file> (default: <input_path>.norm.pnml).\n"
      << "                       All places and transitions have id and name set to p0,p1... and t0,t1...\n"
      << "                       respectively, graphical and tool-specific information is removed.\n"
      << "  -q                   Quiet mode: Suppress detailed invariant output.\n"
      << "  -t <seconds>         Set timeout for computations (default: "
      << DEFAULT_TIMEOUT << "s).\n"
      << "  --noSingleSignRow    Disable single sign row heuristic in invariant computation.\n"
      << "  --noTrivialCull      Disable test for equal or empty columns before algorithm.\n"
      << "  --useCompression     Enable compression for the basis of semiflows.\n"
      << "  --pivot=<strategy>   Set pivot strategy for elimination heuristic:\n"
      << "                       - best: Optimize for best pivot (default).\n"
      << "                       - worst: Use worst pivot (for testing).\n"
      << "                       - first: Use first valid pivot.\n"
      << "  --loopLimit=<n>      Limit elimination loops to <n> iterations (-1 for no limit).\n\n"
      << "Reachability (explicit walk):\n"
      << "  --props=<file>       Load an MCC property XML file (Reachability*.xml).\n"
      << "  --query=<n>          Select property number <n> (0-based, file order); default all.\n"
      << "  --printProps         Print the selected properties (parsed, then normalised) and exit.\n"
      << "  --walkSteps=<n>      Total step budget for the walk (default: until timeout).\n"
      << "  --runLength=<n>      Steps between restarts from the initial marking (default 1000000).\n"
      << "  --seed=<n>           Random seed (default: clock).\n"
      << "  --trace              Print the witness trace as a transition sequence.\n\n"
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

/** Run one random walk toward target; print the MCC verdict line if reached. */
template<typename T>
  bool runWalk (const petri::walk::WalkNet<T> &wnet, const petri::walk::Target<T> &target,
                const std::string &name, const std::string &verdict,
                const petri::walk::WalkBudget &budget, uint64_t seed, bool printTrace, bool quiet)
  {
    petri::walk::RandomStrategy<T> strategy;
    petri::walk::Walker<T> walker (wnet, target, strategy, seed);
    petri::walk::WalkResult res = walker.run (budget);
    const petri::walk::WalkStats &st = res.stats;
    uint64_t ms = st.millis == 0 ? 1 : st.millis;
    std::cout << "Walk " << (res.found ? "found a witness" : "exhausted its budget") << " for " << name
        << " after " << st.steps << " steps, " << st.resets << " resets (" << st.deadEnds
        << " dead ends), " << st.millis << " ms (" << (st.steps / ms) << " steps/ms; "
        << (st.steps ? st.arcVisits / st.steps : 0) << " arc visits and "
        << (st.steps ? st.flips / st.steps : 0) << " flips per step)." << std::endl;
    if (!res.found) return false;
    if (!walker.verify (res.trace)) {
      std::cerr << "Internal error: witness trace does not replay to the goal." << std::endl;
      return false;
    }
    std::cout << "FORMULA " << name << " " << verdict << " TECHNIQUES EXPLICIT RANDOM_WALK" << std::endl;
    if (printTrace) {
      const auto &tnames = wnet.getNet ().getTnames ();
      std::cout << "Witness (" << res.trace.size () << " transitions):";
      for (uint32_t t : res.trace) std::cout << " " << tnames[t];
      std::cout << std::endl;
    } else if (!quiet) {
      std::cout << "Witness length " << res.trace.size () << " ; final marking ";
      walker.currentMarking ().sparse ().print (std::cout);
      std::cout << std::endl;
    }
    return true;
  }

int main_noex (int argc, char *argv[])
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
  std::string exportAsKERSFile;
  bool doNormalize = false;
  std::string normalizePnmlFile;
  std::string loadKERSFile;
  std::string basisKERSFile;
  EliminationHeuristic::PivotStrategy pivotStrategy =
      EliminationHeuristic::PivotStrategy::FindBest;
  ssize_t loopLimit = 500;
  std::string propsFile;
  long query = -1;
  bool printProps = false;
  petri::walk::WalkBudget budget;
  uint64_t seed = static_cast<uint64_t> (std::chrono::steady_clock::now ().time_since_epoch ().count ());
  bool printTrace = false;
  bool doUseQPlusBasis = false;
  bool doUseCompression = false;

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
    } else if (argv[i] == USEQPLUS) {
      doUseQPlusBasis = true;
    } else if (argv[i] == USECOMPRESSION) {
      doUseCompression = true;
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
    } else if (std::string (argv[i]).substr (0, PROPS.size ()) == PROPS) {
      propsFile = std::string (argv[i]).substr (PROPS.size ());
    } else if (std::string (argv[i]).substr (0, QUERY.size ()) == QUERY) {
      query = std::stol (std::string (argv[i]).substr (QUERY.size ()));
    } else if (std::string (argv[i]) == PRINT_PROPS) {
      printProps = true;
    } else if (std::string (argv[i]).substr (0, WALK_STEPS.size ()) == WALK_STEPS) {
      budget.maxSteps = std::stoull (std::string (argv[i]).substr (WALK_STEPS.size ()));
    } else if (std::string (argv[i]).substr (0, RUN_LENGTH.size ()) == RUN_LENGTH) {
      budget.runLength = std::stoull (std::string (argv[i]).substr (RUN_LENGTH.size ()));
    } else if (std::string (argv[i]).substr (0, SEED.size ()) == SEED) {
      seed = std::stoull (std::string (argv[i]).substr (SEED.size ()));
    } else if (std::string (argv[i]) == TRACE) {
      printTrace = true;
    } else if (std::string (argv[i]) == MINFLOWS) {
      minimizeFlows = true;
    } else if (std::string (argv[i]).substr (0, 17) == "--exportAsMatrix=") {
      exportMatrixFile = std::string (argv[i]).substr (17);
    } else if (std::string (argv[i]).substr (0, 15) == "--exportAsKERS=") {
      exportAsKERSFile = std::string (argv[i]).substr (15);
    } else if (std::string (argv[i]).substr (0, 11) == "--loadKERS=") {
      loadKERSFile = std::string (argv[i]).substr (11);
    } else if (std::string (argv[i]).substr (0, 12) == "--basisKERS=") {
      basisKERSFile = std::string (argv[i]).substr (12);
    } else if (std::string (argv[i]).substr (0, 15) == "--normalizePNML") {
      doNormalize = true;  // Set the flag when we see --normalizePNML
      std::string arg = std::string(argv[i]);
      if (arg == NORMALIZE_PNML) {
        normalizePnmlFile = "";  // Empty string means use default later
      } else if (arg.substr(0, 16) == "--normalizePNML=") {
        normalizePnmlFile = arg.substr(16);  // Specific filename provided
      }
    } else {
      std::cout << "[WARNING   ] Option : " << argv[i] << " not recognized"
          << std::endl;
    }
  }

  EliminationHeuristic heur (useSingleSignRow, pivotStrategy, loopLimit,
                             useCulling, minimizeFlows, doUseQPlusBasis, doUseCompression);

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

  // --- KERS matrix-input path (bypass PNML entirely) ---
  if (!loadKERSFile.empty ()) {
    if (!invariants) {
      std::cerr << "Error: --loadKERS requires a flow/semiflow flag (--Pflows, --Psemiflows, --Tflows, --Tsemiflows)." << std::endl;
      return 1;
    }
    if (doUseCompression) {
      std::cerr << "Warning: --useCompression is not supported with --loadKERS (no net structure); flag ignored." << std::endl;
      doUseCompression = false;
    }

    MatrixCol<VAL> M = SparseMatrixIO<VAL>::read (loadKERSFile);
    if (M.getColumnCount () == 0 && M.getRowCount () == 0) {
      return 1;
    }
    bool semi = psemiflows || tsemiflows;
    MatrixCol<VAL> toUse = (tflows || tsemiflows) ? M.transpose () : M;
    EliminationHeuristic heurM (useSingleSignRow, pivotStrategy, loopLimit,
                                useCulling, minimizeFlows, doUseQPlusBasis, false);
    auto time = std::chrono::steady_clock::now ();
    auto [mat, perms] = InvariantMiddle<VAL>::computePInvariants (toUse, semi, timeout, heurM);
    std::string flowKind = (tflows || tsemiflows) ? "T" : "P";
    std::string flowType = semi ? "semi" : "";
    std::cout << "Computed " << mat.getColumnCount () << " " << flowKind << " "
        << flowType << "flows in "
        << std::chrono::duration_cast<std::chrono::milliseconds> (
               std::chrono::steady_clock::now () - time).count () << " ms." << std::endl;
    // basisKERS implies program-to-program mode: skip human-readable output
    bool printHuman = quiet ? false : basisKERSFile.empty ();
    if (printHuman) {
      size_t ndim = M.getRowCount ();
      std::vector<std::string> names;
      names.reserve (ndim);
      for (size_t i = 0; i < ndim; ++i) names.push_back ("x" + std::to_string (i));
      std::vector<VAL> zeros (ndim, 0);
      InvariantMiddle<VAL>::printInvariant (mat, perms, names, zeros);
    }
    if (!basisKERSFile.empty ()) {
      if (!SparseMatrixIO<VAL>::write (mat, basisKERSFile)) return 1;
      std::cout << "Exported basis to " << basisKERSFile << std::endl;
    }
    std::cout << "Total runtime "
        << std::chrono::duration_cast<std::chrono::milliseconds> (
            std::chrono::steady_clock::now () - runtime).count () << " ms." << std::endl;
    return 0;
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

  {
    SparsePetriNet<VAL> *pn = loadXML<VAL> (modelPath);

    // Handle PNML normalization and export only if flag is present
    if (doNormalize) {
      std::string outputFile = normalizePnmlFile;
      if (outputFile.empty()) {
        outputFile = modelPath + ".norm.pnml";
      }
      pn->normalizeNames();
      PNMLExport<VAL>::transform(*pn, outputFile);
    }

    if (draw) {
      std::string title = "Petri Net: " + pn->getName ();
      std::string filename = FlowPrinter<VAL>::drawNet (
          *pn, title, std::set<size_t> (), std::set<size_t> (),
          std::numeric_limits<size_t>::max ());
      std::string targetFile = modelPath + ".dot";
      if (std::rename (filename.c_str (), targetFile.c_str ()) == 0) {
        std::cout << "Renamed output to " << targetFile << std::endl;
      } else {
        std::cerr << "Warning: Could not rename " << filename << " to "
            << targetFile << std::endl;
      }
    }

    // Export matrix if requested (ASCII or KERS binary)
    if (!exportMatrixFile.empty () || !exportAsKERSFile.empty ()) {
      MatrixCol<VAL> sumMatrix = MatrixCol<VAL>::sumProd (-1, pn->getFlowPT (),
                                                          1, pn->getFlowTP ());
      if (!exportMatrixFile.empty ())
        MatrixExporter<VAL>::exportMatrix (sumMatrix, exportMatrixFile);
      if (!exportAsKERSFile.empty ())
        SparseMatrixIO<VAL>::write (sumMatrix, exportAsKERSFile);
    }

    if (!propsFile.empty ()) {
      std::vector<petri::expr::Property> props = petri::mcc::loadProperties (propsFile, *pn);
      if (query >= 0) {
        if (static_cast<size_t> (query) >= props.size ()) {
          std::cerr << "Error: --query=" << query << " but the file holds " << props.size () << " properties." << std::endl;
          delete pn;
          return 1;
        }
        props = { props[static_cast<size_t> (query)] };
      }
      if (printProps) {
        for (const auto &prop : props) {
          prop.print (std::cout, &pn->getPnames ());
          std::cout << "\n  goal (" << prop.goal ().size () << " nodes) : ";
          prop.goal ().print (std::cout, &pn->getPnames ());
          std::cout << std::endl;
        }
        delete pn;
        return 0;
      }
      budget.timeoutMillis = static_cast<uint64_t> (timeout) * 1000;
      petri::walk::WalkNet<VAL> wnet (*pn);
      for (const auto &prop : props) {
        if (prop.kind == petri::expr::PropertyKind::Unsupported) {
          std::cout << "Skipping " << prop.name << " : " << prop.comment << std::endl;
          continue;
        }
        petri::walk::Target<VAL> target = prop.kind == petri::expr::PropertyKind::Deadlock
            ? petri::walk::Target<VAL>::deadlockTarget ()
            : petri::walk::Target<VAL> (prop.goal ());
        if (!target.isDeadlock () && target.expression ().kind == petri::expr::Expression::Kind::False) {
          std::cout << "FORMULA " << prop.name << " " << (prop.kind == petri::expr::PropertyKind::Invariant ? "TRUE" : "FALSE")
              << " TECHNIQUES TOPOLOGICAL TRIVIAL" << std::endl;
          continue;
        }
        runWalk<VAL> (wnet, target, prop.name, prop.verdictIfReached (), budget, seed, printTrace, quiet);
      }
    }

    if (findDeadlock) {
      budget.timeoutMillis = static_cast<uint64_t> (timeout) * 1000;
      petri::walk::WalkNet<VAL> wnet (*pn);
      petri::walk::Target<VAL> target = petri::walk::Target<VAL>::deadlockTarget ();
      runWalk<VAL> (wnet, target, "ReachabilityDeadlock", "TRUE", budget, seed, printTrace, quiet);
    }

    if (invariants) {
      vector<int> repr;
      if (pflows || psemiflows) {
        auto time = std::chrono::steady_clock::now ();
        MatrixCol<VAL> sumMatrix = MatrixCol<VAL>::sumProd (-1,
                                                            pn->getFlowPT (), 1,
                                                            pn->getFlowTP ());
        auto [mat,perms] = InvariantMiddle<VAL>::computePInvariants (sumMatrix,
                                                               psemiflows,
                                                               timeout, heur);
        std::cout << "Computed " << mat.getColumnCount () << " P "
            << (psemiflows ? "semi" : "") << "flows " ;
        if (!perms.empty()) {
          std::cout << "with " << perms.size() << " permutations.\n";
          std::cout << "Total decompressed invariants :" << InvariantMiddle<VAL>::countInvariant(mat,perms) << std::endl;
        }
        std::cout << " in " << std::chrono::duration_cast<std::chrono::milliseconds> (
                std::chrono::steady_clock::now () - time).count () << " ms."
            << std::endl;
        if (!quiet && basisKERSFile.empty ()) {
          InvariantMiddle<VAL>::printInvariant (mat, perms, pn->getPnames (),
                                                (*pn).getMarks ());
        }
        if (!basisKERSFile.empty ()) {
          SparseMatrixIO<VAL>::write (mat, basisKERSFile);
          std::cout << "Exported basis to " << basisKERSFile << std::endl;
        }
      }
      if (tflows || tsemiflows) {
        auto time = std::chrono::steady_clock::now ();
        MatrixCol<VAL> sumMatrix =
            MatrixCol<VAL>::sumProd (-1, pn->getFlowPT (), 1, pn->getFlowTP ()).transpose ();
        auto [mat,perms] = InvariantMiddle<VAL>::computePInvariants (sumMatrix,
                                                               tsemiflows,
                                                               timeout, heur);
        std::cout << "Computed " << mat.getColumnCount () << " T "
            << (tsemiflows ? "semi" : "") << "flows ";
        if (!perms.empty()) {
          std::cout << "with " << perms.size() << " permutations.\n";
          std::cout << "Total decompressed invariants :" << InvariantMiddle<VAL>::countInvariant(mat,perms) << std::endl;
        }
        std::cout << " in " << std::chrono::duration_cast<std::chrono::milliseconds> (
                std::chrono::steady_clock::now () - time).count () << " ms."
            << std::endl;
        if (!quiet && basisKERSFile.empty ()) {
          std::vector<VAL> emptyVector (pn->getTnames ().size (), 0);
          InvariantMiddle<VAL>::printInvariant (mat, perms, pn->getTnames (),
                                                emptyVector);
        }
        if (!basisKERSFile.empty ()) {
          SparseMatrixIO<VAL>::write (mat, basisKERSFile);
          std::cout << "Exported basis to " << basisKERSFile << std::endl;
        }
      }
    }

    delete pn;
  }

  std::cout << "Total runtime "
      << std::chrono::duration_cast<std::chrono::milliseconds> (
          std::chrono::steady_clock::now () - runtime).count () << " ms."
      << std::endl;

  return 0;
}

int main (int argc, char *argv[])
{
  // Reserve 16K that can be freed to allow an OOM error message to be printed
  char *emergencyMemory = new char[16384];
  try {
    return main_noex (argc, argv);
  } catch (const char *ex) {
    std::cerr << "An unexpected exception occurred : " << ex << std::endl;
    return 1;
  } catch (std::string &ex) {
    std::cerr << "An unexpected exception occurred : " << ex << std::endl;
    return 1;
  } catch (std::bad_alloc &) {
    delete[] emergencyMemory;
    std::cerr << "Out of memory error!" << std::endl;
    return 1;
  }
}
