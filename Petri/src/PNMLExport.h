/*
 * PNMLExport.h
 *
 * Header file for exporting SparsePetriNet to PNML format
 */

#ifndef PNMLEXPORT_H_
#define PNMLEXPORT_H_

#include <fstream>
#include <string>
#include <chrono>
#include <iostream>
#include "SparsePetriNet.h"

template<typename T>
  class PNMLExport
  {
  public:
    static void transform (const SparsePetriNet<T> &sr, const std::string &path)
    {
      // Start timing
      auto start = std::chrono::high_resolution_clock::now ();

      // Open output file
      std::ofstream pw (path);
      if (!pw.is_open ()) {
        throw std::runtime_error ("Unable to open file: " + path);
      }

      // Write PNML header
      pw << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
      pw << "<pnml xmlns=\"http://www.pnml.org/version-2009/grammar/pnml\">\n";
      pw << "<net id=\"" << sr.getName ()
          << "\" type=\"http://www.pnml.org/version-2009/grammar/ptnet\">\n";
      pw << "<page id=\"page0\"><name><text>DefaultPage</text></name>\n";

      // Export places
      int nbplaces = 0;
      for (size_t p = 0; p < sr.getPnames ().size (); p++) {
        pw << "<place id=\"p" << p << "\">";
        pw << "<name><text>" << sr.getPnames ()[p] << "</text></name>";
        T mark = sr.getMarks ()[p];
        if (mark != 0) {
          pw << "<initialMarking><text>" << mark << "</text></initialMarking>";
        }
        pw << "</place>\n";
        nbplaces++;
      }

      // Export transitions
      int nbtrans = 0;
      for (size_t t = 0; t < sr.getTnames ().size (); t++) {
        pw << "<transition id=\"t" << t << "\">";
        pw << "<name><text>" << sr.getTnames ()[t] << "</text></name>";
        pw << "</transition>\n";
        nbtrans++;
      }

      // Export arcs
      int arcid = 0;
      // PT arcs
      for (size_t t = 0; t < sr.getTnames ().size (); t++) {
        const SparseArray<T> &pt = sr.getFlowPT ().getColumn (t);
        for (size_t i = 0; i < pt.size (); i++) {
          size_t p = pt.keyAt (i);
          T val = pt.valueAt (i);
          pw << "<arc id=\"a" << (arcid++) << "\" source=\"p" << p
              << "\" target=\"t" << t << "\">";
          if (val != 1) {
            pw << "<inscription><text>" << val << "</text></inscription>";
          }
          pw << "</arc>\n";
        }
      }

      // TP arcs
      for (size_t t = 0; t < sr.getTnames ().size (); t++) {
        const SparseArray<T> &tp = sr.getFlowTP ().getColumn (t);
        for (size_t i = 0; i < tp.size (); i++) {
          size_t p = tp.keyAt (i);
          T val = tp.valueAt (i);
          pw << "<arc id=\"a" << (arcid++) << "\" source=\"t" << t
              << "\" target=\"p" << p << "\">";
          if (val != 1) {
            pw << "<inscription><text>" << val << "</text></inscription>";
          }
          pw << "</arc>\n";
        }
      }

      // Close tags
      pw << "</page>\n";
      pw << "<name><text>" << sr.getName () << "</text></name></net>\n";
      pw << "</pnml>\n";

      // Flush and close file
      pw.flush ();
      pw.close ();

      // Calculate and log timing
      auto end = std::chrono::high_resolution_clock::now ();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds> (
          end - start).count ();

      std::cout << "Export to PNML in file " << path << " of net with "
          << nbplaces << " places, " << nbtrans << " transitions and " << arcid
          << " arcs took " << duration << " ms." << std::endl;
    }
  };

#endif /* PNMLEXPORT_H_ */
