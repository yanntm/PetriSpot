/*
 * SparsePetriNet.h
 *
 *  Created on: Oct 13, 2020
 *      Author: ythierry
 */

#ifndef SPARSEPETRINET_H_
#define SPARSEPETRINET_H_

#include "MatrixCol.h"

template<typename T>
  class SparsePetriNet
  {
    std::string name;
    std::vector<T> marks;
    MatrixCol<T> flowPT;
    MatrixCol<T> flowTP;
    std::vector<std::string> tnames;
    std::vector<std::string> pnames;
    T maxArcValue;
    static const int DEBUG = 0;
  public:
    SparsePetriNet ()
        : name ("Petri"), maxArcValue (0)
    {
    }

    const std::string& getName () const
    {
      return this->name;
    }

    void setName (const std::string &n)
    {
      this->name = n;
    }

    size_t addTransition (const std::string &tname)
    {
      flowPT.appendColumn (SparseArray<T> ());
      flowTP.appendColumn (SparseArray<T> ());
      tnames.emplace_back (tname);
      return tnames.size () - 1;
    }

    size_t addPlace (const std::string &pname, T init)
    {
      flowPT.addRow ();
      flowTP.addRow ();
      pnames.emplace_back (pname);
      marks.push_back (init);
      return pnames.size () - 1;
    }

    void addPreArc (int p, int t, T val)
    {
      flowPT.getColumn (t).put (p, val);
      maxArcValue = std::max (maxArcValue, val);
    }

    void addPostArc (int p, int t, T val)
    {
      flowTP.getColumn (t).put (p, val);
      maxArcValue = std::max (maxArcValue, val);
    }

    size_t getTransitionCount () const
    {
      return tnames.size ();
    }

    size_t getPlaceCount () const
    {
      return pnames.size ();
    }

    size_t getArcCount () const
    {
      size_t sum = 0;
      for (size_t i = 0; i < flowTP.getColumnCount (); i++) {
        sum += flowTP.getColumn (i).size ();
      }
      for (size_t i = 0; i < flowTP.getColumnCount (); i++) {
        sum += flowPT.getColumn (i).size ();
      }
      return sum;
    }

    void setMarking (size_t pid, T val)
    {
      marks[pid] = val;
    }

    const std::vector<std::string>& getTnames () const
    {
      return tnames;
    }

    const std::vector<std::string>& getPnames () const
    {
      return pnames;
    }

    void normalizeNames () {
      for (auto &name : tnames) {
        normalizeName (name);
      }
      for (auto &name : pnames) {
        normalizeName (name);
      }
    }

    static void normalizeName (std::string &name) {
      std::replace (name.begin (), name.end (), ' ', '_');
      std::replace (name.begin (), name.end (), '-', '_');
      std::replace (name.begin (), name.end (), '/', '_');
      std::replace (name.begin (), name.end (), '*', 'x');
      std::replace (name.begin (), name.end (), '=', '_');
    }

    MatrixCol<T>& getFlowPT ()
    {
      return flowPT;
    }
    MatrixCol<T>& getFlowTP ()
    {
      return flowTP;
    }
    const MatrixCol<T>& getFlowPT () const
    {
      return flowPT;
    }
    const MatrixCol<T>& getFlowTP () const
    {
      return flowTP;
    }

    int getMaxArcValue ()
    {
      return maxArcValue;
    }

    const std::vector<T>& getMarks () const
    {
      return marks;
    }
  };

#endif /* SPARSEPETRINET_H_ */
