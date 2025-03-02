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

namespace petri
{
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

    static void printInvariant (const MatrixCol<T> &invariants,
                                const std::vector<std::string> &pnames,
                                const std::vector<T> &initial,
                                std::ostream &out)
    {
      for (const auto &rv : invariants.getColumns ()) {
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
      out << "Total of " << invariants.getColumnCount () << " invariants."
          << std::endl;
    }

    static void printInvariant (const MatrixCol<T> &invariants,
                                const std::vector<std::string> &pnames,
                                const std::vector<T> &initial)
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
        size_t k = inv.keyAt (i);
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

    static MatrixCol<T> computePInvariants (const MatrixCol<T> &pn)
    {
      return computePInvariants (pn, false, 120);
    }

    static MatrixCol<T> computePInvariants (
        MatrixCol<T> &pn, bool onlyPositive, int timeout,
        const EliminationHeuristic &heuristic = EliminationHeuristic ())
    {
      std::promise<MatrixCol<T> > promise;
      auto future = promise.get_future ();

      std::thread thread ( [&] () {
        promise.set_value(computePInvariants(pn, onlyPositive, heuristic));
      });

      auto status = future.wait_for (std::chrono::seconds (timeout));
      if (status == std::future_status::ready) {
        thread.join ();
        return future.get ();
      } else {
        if (thread.joinable ()) {
          thread.detach ();
        }
        writeToLog (
            "Time Limit : exiting after meeting timeout of "
                + std::to_string (timeout) + " seconds.");
        return MatrixCol<T> ();
      }
    }

  public:

    static MatrixCol<T> computePInvariants (
        MatrixCol<T> &pn, bool onlyPositive,
        const EliminationHeuristic &heuristic = EliminationHeuristic ())
    {
      auto startTime = std::chrono::steady_clock::now ();
      try {
        auto tpn = pn.transpose ();
        auto mat = petri::InvariantCalculator<T>::calcInvariantsPIPE (
            tpn, onlyPositive, heuristic);
        std::string logMessage = "Computed "
            + std::to_string (mat.getColumnCount ()) + " invariants in "
            + std::to_string (
                std::chrono::duration_cast<std::chrono::milliseconds> (
                    std::chrono::steady_clock::now () - startTime).count ())
            + " ms";
        writeToLog (logMessage);
        return mat;
      } catch (std::overflow_error &e) {
        std::cerr << e.what () << std::endl;
        std::string logMessage = "Invariants computation overflowed in "
            + std::to_string (
                std::chrono::duration_cast<std::chrono::milliseconds> (
                    std::chrono::steady_clock::now () - startTime).count ())
            + " ms";
        writeToLog (logMessage);
        return MatrixCol<T> ();
      }
    }

  };

} /* namespace petri */

#endif /* INVARIANTMIDDLE_H_ */
