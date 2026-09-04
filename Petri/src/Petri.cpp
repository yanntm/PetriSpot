/*
 * Petri.cpp
 *
 * Entry point of the petri32/64/128 binaries: parse the command line
 * (cli/Options.h), load the input, dispatch to the invariant and walk drivers.
 */
#include <chrono>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <string>

#include "cli/InvariantDriver.h"
#include "cli/Options.h"
#include "cli/WalkDriver.h"
#include "cli11/CLI11.hpp"
#include "core/Log.h"
#include "core/SparsePetriNet.h"
#include "io/FlowPrinter.h"
#include "io/MatrixExporter.h"
#include "io/PNETIO.h"
#include "io/PNMLExport.h"
#include "io/SparseMatrixIO.h"
#include "parse/PTNetLoader.h"
#include "walk/NetStats.h"
#include "walk/WalkNet.h"

#ifndef VAL
// default to 64 bit
#define VAL long
#endif

using petri::cli::Options;

/** Exports and inspections of a loaded net requested by the options. */
static void exportNet (const Options &o, SparsePetriNet<VAL> &pn)
{
  if (o.normalizePNML) {
    std::string base = o.modelPath.empty () ? o.netFile : o.modelPath;
    std::string out = o.normalizePNMLFile.empty () ? base + ".norm.pnml" : o.normalizePNMLFile;
    pn.normalizeNames ();
    PNMLExport<VAL>::transform (pn, out);
  }
  if (o.draw) {
    std::string filename = FlowPrinter<VAL>::drawNet (pn, "Petri Net: " + pn.getName (), std::set<size_t> (),
                                                      std::set<size_t> (), std::numeric_limits<size_t>::max ());
    std::string target = (o.modelPath.empty () ? o.netFile : o.modelPath) + ".dot";
    if (std::rename (filename.c_str (), target.c_str ()) == 0) {
      std::cout << "Renamed output to " << target << std::endl;
    } else {
      std::cerr << "Warning: Could not rename " << filename << " to " << target << std::endl;
    }
  }
  if (!o.exportMatrixFile.empty () || !o.exportAsKERSFile.empty ()) {
    MatrixCol<VAL> sumMatrix = MatrixCol<VAL>::sumProd (-1, pn.getFlowPT (), 1, pn.getFlowTP ());
    if (!o.exportMatrixFile.empty ()) MatrixExporter<VAL>::exportMatrix (sumMatrix, o.exportMatrixFile);
    if (!o.exportAsKERSFile.empty ()) SparseMatrixIO<VAL>::write (sumMatrix, o.exportAsKERSFile);
  }
  if (!o.exportNetFile.empty () && PNETIO<VAL>::write (pn, o.exportNetFile)) {
    std::cout << "Exported net to " << o.exportNetFile << std::endl;
  }
  if (o.netStats) {
    petri::walk::WalkNet<VAL> wnet (pn);
    petri::walk::printNetStats (std::cout, wnet);
  }
}

static int main_noex (int argc, char *argv[])
{
  std::string logMessage = "Running PetriSpot with arguments : [";
  for (int i = 1; i < argc; ++i) logMessage += std::string (argv[i]) + (i + 1 < argc ? ", " : "");
  petri::writeToLog (logMessage + (argc == 1 ? "Empty" : "") + "]");
  auto runtime = std::chrono::steady_clock::now ();

  Options o;
  CLI::App app ("PetriSpot: invariants and explicit reachability walks for P/T Petri nets.");
  petri::cli::addOptions (app, o);
  if (argc == 1) {
    std::cerr << app.help () << std::endl;
    return 1;
  }
  try {
    app.parse (argc, argv);
  } catch (const CLI::ParseError &e) {
    return app.exit (e);
  }
  if (o.pflows && o.psemiflows) {
    std::cout << "Cannot compute P flows and P semi-flows at the same time." << std::endl;
    return 1;
  }
  if (o.tflows && o.tsemiflows) {
    std::cout << "Cannot compute T flows and T semi-flows at the same time." << std::endl;
    return 1;
  }

  int status = 0;
  if (!o.loadKERSFile.empty ()) {
    status = petri::cli::runInvariantsFromKERS<VAL> (o);
  } else {
    if (o.modelPath.empty () && o.netFile.empty ()) {
      std::cerr << "Error: no model file specified (-i or --net)." << std::endl;
      return 1;
    }
    std::unique_ptr<SparsePetriNet<VAL>> pn;
    if (!o.netFile.empty ()) {
      auto time = std::chrono::steady_clock::now ();
      pn.reset (PNETIO<VAL>::read (o.netFile));
      petri::writeToLog ("Loaded PNET net with " + std::to_string (pn->getPlaceCount ()) + " places, "
          + std::to_string (pn->getTransitionCount ()) + " transitions and " + std::to_string (pn->getArcCount ())
          + " arcs in " + std::to_string (petri::cli::millisSince (time)) + " ms.");
    } else {
      if (!std::ifstream (o.modelPath).good ()) {
        std::cerr << "Error: file not found: " << o.modelPath << std::endl;
        return 1;
      }
      pn.reset (loadXML<VAL> (o.modelPath));
    }
    exportNet (o, *pn);
    if (!o.propsFile.empty ()) petri::cli::runProperties (o, *pn);
    if (o.findDeadlock) petri::cli::runDeadlock (o, *pn);
    if (o.invariants ()) petri::cli::runInvariants (o, *pn);
  }
  std::cout << "Total runtime " << petri::cli::millisSince (runtime) << " ms." << std::endl;
  return status;
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
