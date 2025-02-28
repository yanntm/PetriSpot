#ifndef INVARIANTMIDDLE_H_
#define INVARIANTMIDDLE_H_

#include "Arithmetic.hpp"
#include "SparsePetriNet.h"
#include "InvariantCalculator.h"
#include "Heuristic.h"
#include <chrono>
#include <thread>
#include <future>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <unordered_map>
#include <unordered_set>


namespace petri {
/**
 * A front-end for functionality computing invariants. Underlying code is
 * adapted from CvO-Theory group's APT : https://github.com/CvO-Theory/apt 
 * See also uniol.apt package and classes.
 * 
 * @author ythierry
 *
 */
template<typename T>
  class InvariantMiddle
  {
  private:
    static inline MatrixCol<T> last;
    static inline std::mutex lock;
    static inline std::unordered_set<SparseArray<T>> lastInv;
  public:
    static const int DEBUG = 0;

    /**
     * Guaranteed polynomial runtime, returns flows (with positive AND negative
     * coefficients)
     *
     * @param pn representing the Petri net approximation
     * @return a set of invariants, i.e. coeffs for each variable such that the sum
     *         is constant in all markings/states.
     */

    /**
     * Worst case exponential (time and memory), returns semi-flows (with positive
     * coefficients only) which are reputed easier to interpret.
     *
     * @param pn representing the Petri net approximation
     * @return a set of invariants, i.e. coeffs for each variable such that the sum
     *         is constant in all markings/states.
     */


    static void printInvariant (
        const std::unordered_set<SparseArray<T>> &invariants,
        const std::vector<std::string> &pnames, const std::vector<T> &initial,
        std::ostream &out)
    {
      for (const auto &rv : invariants) {
        std::stringstream sb;
        try {
          auto sum = printEquation (rv, initial, pnames, sb);
          out << "inv : " << sb.str () << " = " << sum << "\n";
        } catch (std::overflow_error &e) {
          std::cerr
              << "Overflow of 'int' when computing constant for invariant."
              << std::endl;
        }
      }
      out << "Total of " << invariants.size () << " invariants." << std::endl;
    }

    static void printInvariant (
        const std::unordered_set<SparseArray<T>> &invariants,
        const std::vector<std::string> &pnames, const std::vector<T> &initial)
    {
      printInvariant (invariants, pnames, initial, std::cout);
    }

    static T printEquation (const SparseArray<T> &inv,
                              const std::vector<T> &initial,
                              const std::vector<std::string> &pnames,
                              std::stringstream &sb)
    {
      bool first = true;
      T sum = 0;
      for (size_t i = 0; i < inv.size (); i++) {
        unsigned int k = inv.keyAt (i);
        T v = inv.valueAt (i);
        if (v != 0) {
          if (!initial.empty ()) {
             sum = petri::addExact (sum, (petri::multiplyExact (v, initial[k])));
          }
          if (!first) {
            if (v < 0) {
              sb << " - ";
              v = -v;
            } else {
              sb << " + ";
            }
          } else {
            if (v < 0) {
              sb << "-";
              v = -v;
            }
            first = false;
          }
          if (v != 1) {
            sb << v << "*" << pnames[k];
          } else {
            sb << pnames[k];
          }
        }
      }
      return sum;
    }

  public:
    static void writeToLog (const std::string &message)
    {
      auto now = std::chrono::system_clock::now ();
      std::time_t now_c = std::chrono::system_clock::to_time_t (now);
      struct std::tm *parts = std::localtime (&now_c);

      std::cout << "[" << (parts->tm_year + 1900) << "-" << std::setw (2)
          << std::setfill ('0') << (parts->tm_mon + 1) << "-" << std::setw (2)
          << std::setfill ('0') << parts->tm_mday << " " << std::setw (2)
          << std::setfill ('0') << parts->tm_hour << ":" << std::setw (2)
          << std::setfill ('0') << parts->tm_min << ":" << std::setw (2)
          << std::setfill ('0') << parts->tm_sec << "] " << "[INFO   ] "
          << message << std::endl;
    }

    static std::unordered_set<SparseArray<T>> computePInvariants (
        const MatrixCol<T> &pn)
    {
      return computePInvariants (pn, false, 120);
    }

    static std::unordered_set<SparseArray<T>> computePInvariants (
        const MatrixCol<T> &pn, bool onlyPositive, int timeout, const EliminationHeuristic &heuristic=EliminationHeuristic())
    {
      std::promise < std::unordered_set<SparseArray<T>> > promise;
      auto future = promise.get_future ();

      std::thread thread (
          [&] () {
            std::unordered_set<SparseArray<T>> result = computePInvariants(pn, onlyPositive, heuristic);
            promise.set_value(result);
          });

      auto status = future.wait_for (std::chrono::seconds (timeout));
      if (status == std::future_status::ready) {
        thread.join ();
        return future.get ();
      } else {
        if (thread.joinable ()) {
          thread.detach ();
        }
        writeToLog ("Time Limit : exiting after meeting timeout of "+std::to_string(timeout)+" seconds.");
        return std::unordered_set<SparseArray<T>> ();
      }
    }

  private:
    static void cache (const MatrixCol<T> &pn,
                       const std::unordered_set<SparseArray<T>> &inv)
    {
      std::lock_guard < std::mutex > guard (lock);
      last = pn;
      lastInv = inv;
    }

    static std::unordered_set<SparseArray<T>> checkCache (
        const MatrixCol<T> &pn)
    {
      std::lock_guard < std::mutex > guard (lock);
      if (pn == last) {
        writeToLog ("Invariant cache hit.");
        return lastInv;
      } else {
        return std::unordered_set<SparseArray<T>> ();
      }
    }

  public:
    static std::unordered_set<SparseArray<T>> computeTinvariants (
        const SparsePetriNet<T> &sr, const MatrixCol<T> &sumMatrix,
        const std::vector<int> &repr, bool onlyPositive)
    {

      std::unordered_map<int, std::vector<int>> repSet = computeMap (repr);
      std::unordered_set<SparseArray<T>> invarT = computePInvariants (
          sumMatrix.transpose (), onlyPositive);

      if (DEBUG >= 1 && !invarT.empty ()) {
        std::vector < T > empty (sr.getTransitionCount (), 0);
        printInvariant (invarT, sr.getTnames (), empty);
      }
      // so we have T invariants, using the reduced flow matrix
      std::unordered_set<SparseArray<T>> reindexT;
      // reinterpret over the original indexes of transitions
      for (const auto &inv : invarT) {
        std::vector<SparseArray<T>> toadd =
          { SparseArray<T> () };
        for (size_t i = 0; i < inv.size (); ++i) {
          int t = inv.keyAt (i);
          int val = inv.valueAt (i);
          const std::vector<int> &images = repSet[t];
          if (images.size () > 1) {
            for (const auto &img : images) {
              std::vector<SparseArray<T>> toadd2;
              for (const auto &b : toadd) {
                SparseArray<T> mod = b;
                mod.put (img, val);
                toadd2.push_back (mod);
              }
              toadd = std::move (toadd2);
            }
          } else {
            for (auto &b : toadd) {
              b.put (images[0], val);
            }
          }
        }
        for (const auto &b : toadd) {
          reindexT.insert (b);
        }
      }

      if (DEBUG >= 2 && !invarT.empty ()) {
        printInvariant (reindexT, sr.getTnames (),
          { });
      }

      return invarT;
    }

    static std::unordered_set<SparseArray<T>> computeTinvariants (
        const SparsePetriNet<T> &sr, MatrixCol<T> &sumMatrix,
        const std::vector<int> &repr)
    {
      return computeTinvariants (sr, sumMatrix, repr, true);
    }

    static std::unordered_set<SparseArray<T>> computePInvariants (
        const MatrixCol<T> &pn, bool onlyPositive,  const EliminationHeuristic &heuristic=EliminationHeuristic())
    {
      std::unordered_set<SparseArray<T>> invar = checkCache (pn);
      if (!invar.empty ()) {
        return invar;
      }

      auto startTime = std::chrono::steady_clock::now ();
      try {
        invar = petri::InvariantCalculator<T>::calcInvariantsPIPE (pn.transpose (),
                                                            onlyPositive, heuristic);
        cache (pn, invar);
        std::string logMessage = "Computed " + std::to_string (invar.size ())
            + " invariants in "
            + std::to_string (
                std::chrono::duration_cast < std::chrono::milliseconds
                    > (std::chrono::steady_clock::now () - startTime).count ())
            + " ms";
        writeToLog (logMessage);
      } catch (std::overflow_error &e) {
        std::cerr << e.what () << std::endl;
        invar.clear ();
        std::string logMessage = "Invariants computation overflowed in "
            + std::to_string (
                std::chrono::duration_cast < std::chrono::milliseconds
                    > (std::chrono::steady_clock::now () - startTime).count ())
            + " ms";
        writeToLog (logMessage);
      }
      return invar;
    }



    /**
     * Computes a combined flow matrix, stored with column = transition, while
     * removing any duplicates (e.g. due to test arcs or plain redundancy). Updates
     * tnames that is supposed to initially be empty to set the names of the
     * transitions that were kept. This is so we can reinterpret appropriately the
     * Parikh vectors f so desired.
     *
     * @param sr             our Petri net
     * @param representative the mapping from original transition index to their new
     *                       representative (many to one/surjection)
     * @return a (reduced, less columns than usual) flow matrix
     */
    static MatrixCol<T> computeReducedFlow (const SparsePetriNet<T> &sr,
                                            std::vector<int> &representative)
    {
      MatrixCol<T> sumMatrix (sr.getPlaceCount (), 0);
      sumMatrix.reserveColumns (sr.getTransitionCount ());
      {
        typedef std::unordered_map<const SparseArray<T>*, int,
            std::hash<SparseArray<T>*>, std::equal_to<SparseArray<T>*>> map_t;
        map_t seen;
        int discarded = 0;
        int curr = 0;

        for (size_t i = 0; i < sr.getTransitionCount (); i++) {
          // effects of t
          SparseArray<T> combined = SparseArray<T>::sumProd (
              -1, sr.getFlowPT ().getColumn (i), 1,
              sr.getFlowTP ().getColumn (i));
          // have we seen this effect ?
          typename map_t::iterator it = seen.find (&combined);
          if (it == seen.end ()) {
            // a new effect
            sumMatrix.appendColumn (combined);
            seen.insert (
                std::make_pair (
                    &sumMatrix.getColumn (sumMatrix.getColumnCount () - 1),
                    curr));
            // this transition is its own representative
            representative.push_back (curr);
            curr++;
          } else {
            // this transition is represented by the element at index :
            representative.push_back (it->second);
            discarded++;
          }
        }
        if (discarded > 0) {
          std::string logMessage = "Flow matrix only has "
              + std::to_string (sumMatrix.getColumnCount ())
              + " transitions (discarded " + std::to_string (discarded)
              + " similar events)";
          writeToLog (logMessage);
        }
      }
      return sumMatrix;
    }

    static SparseArray<T> transformParikh (
        const SparseArray<T> &parikhori,
        const std::unordered_map<int, std::vector<int>> &repr)
    {
      SparseArray<T> parikh;
      for (int i = 0, e = parikhori.size (); i < e; i++) {
        int t = parikhori.keyAt (i);
        int k = parikhori.valueAt (i);
        auto it = repr.find (t);
        if (it != repr.end ()) {
          for (int tr : it->second) {
            parikh.put (tr, k);
          }
        }
      }
      return parikh;
    }

    static std::unordered_map<int, std::vector<int>> computeMap (
        const std::vector<int> &repr)
    {
      std::unordered_map<int, std::vector<int>> repSet;
      for (size_t i = 0; i < repr.size (); ++i) {
        int t = i;
        repSet[repr[t]].push_back (t);
      }
      return repSet;
    }

  };

} /* namespace petri */

#endif /* INVARIANTMIDDLE_H_ */
