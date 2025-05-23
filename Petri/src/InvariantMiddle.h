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

    using Permutations = typename InvariantCalculator<T>::Permutations;
    using Invariants = std::pair<MatrixCol<T>, Permutations>;

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
                                const Permutations &permutations,
                                const std::vector<std::string> &pnames,
                                const std::vector<T> &initial,
                                std::ostream &out)
    {
      std::vector<std::string> moreNames;
      std::vector<T> moreInitials;

      size_t firstPerm = pnames.size ();

      std::unordered_map<size_t,__uint128_t> repSizes;
      if (!permutations.empty ()) {
        out << "Permutations : " << std::endl;

        // generate names
        moreNames.reserve (pnames.size () * permutations.size ());
        std::copy (pnames.begin (), pnames.end (), std::back_inserter (moreNames));
        moreInitials.reserve (pnames.size () * permutations.size ());
        std::copy (initial.begin (), initial.end (), std::back_inserter (moreInitials));

        // Now print the permutations
        for (size_t i = 0; i < permutations.size (); i++) {
          // format is :
          // r0 : [ p0 + 2*p1 = 2 ] , [ p3 = 1 ] ; (2)
          auto name = "r" + std::to_string (i);
          moreNames.push_back(name);
          moreInitials.push_back(0);
          out << name << " : ";

          const auto & perm = permutations[i];

          __uint128_t represents = 0;

          bool first = true;
          for (const auto &term : perm.elements) {

            if (!first) {
              out << ", ";
            } else {
              first = false;
            }
            out << "[ " ;
            std::stringstream sb;
            auto sum = printEquation (term, moreInitials, moreNames, sb);
            out << sb.str () << " = " << sum << " ]";
            __uint128_t termrepresents = 1;
            // we need to add constants to right hand side
            for (size_t i = term.size (); i-->0 ;) {
              size_t key = term.keyAt (i);
              if (key < firstPerm) {
                break;
              } else {
                termrepresents = petri::multiplyExact(termrepresents,repSizes[key]);
              }
            }
            represents += termrepresents ;
          }
          repSizes[perm.index] = represents;

          out << " ; (" << represents << ") \n";
        }
      }
      // Declare and initialize names and initials conditionally
      const std::vector<std::string> &names = permutations.empty() ? pnames : moreNames;
      const std::vector<T> &initials = permutations.empty() ? initial : moreInitials;

      out << "Invariants : \n" ;
      __uint128_t decompressed = 0;
      for (const auto &rv : invariants.getColumns ()) {
        std::stringstream sb;
        try {
          auto sum = printEquation (rv, initials, names, sb);
          out << "inv : " << sb.str () << " = " << sum ;
          __uint128_t represents = 1;
          if (!permutations.empty ()) {
            // we need to add constants to right hand side
            for (size_t i = rv.size (); i-->0 ;) {
              size_t key = rv.keyAt (i);
              if (key < firstPerm) {
                break;
              } else {
                T val = rv.valueAt (i);
                out << (val < 0 ? " - " : " + ") << (std::abs(val)==1?"":std::to_string(val)+"*") << "kr" << (key-firstPerm);
                represents = petri::multiplyExact(represents,repSizes[key]);
              }
            }
            out << " (" << represents << ")";
            out << "\n";
          } else {
            out << "\n";
          }
          decompressed = petri::addExact(decompressed,represents);
        } catch (std::overflow_error &e) {
          std::cerr
              << "Overflow of 'int' when computing constant for invariant."
              << std::endl;
        }
      }
      out << "Total of " << invariants.getColumnCount () << (permutations.empty()?"":" compressed") << " invariants."
          << std::endl;
      if (!permutations.empty ()) {
        out << "Total of " << decompressed << " decompressed invariants."
            << std::endl;
      }
    }

    static __uint128_t countInvariant (const MatrixCol<T> &invariants,
                                const Permutations &permutations)
    {
      if (permutations.empty ()) {
        return invariants.getColumnCount ();
      }

      size_t firstPerm = invariants.getRowCount () - permutations.size ();

      std::unordered_map<size_t,__uint128_t> repSizes;
      if (!permutations.empty ()) {
        // Now print the permutations
        for (size_t i = 0; i < permutations.size (); i++) {
          // format is :
          // r0 : [ p0 + 2*p1 = 2 ] , [ p3 = 1 ] ; (2)
          const auto & perm = permutations[i];

          __uint128_t represents = 0;

          for (const auto &term : perm.elements) {
            __uint128_t termrepresents = 1;
            // we need to add constants to right hand side
            for (size_t i = term.size (); i-->0 ;) {
              size_t key = term.keyAt (i);
              if (key < firstPerm) {
                break;
              } else {
                termrepresents = petri::multiplyExact(termrepresents,repSizes[key]);
              }
            }
            represents += termrepresents ;
          }
          repSizes[perm.index] = represents;
        }
      }
      __uint128_t decompressed = 0;
      for (const auto &rv : invariants.getColumns ()) {
        try {
          __uint128_t represents = 1;
          if (!permutations.empty ()) {
            // we need to add constants to right hand side
            for (size_t i = rv.size (); i-->0 ;) {
              size_t key = rv.keyAt (i);
              if (key < firstPerm) {
                break;
              } else {
                represents = petri::multiplyExact(represents,repSizes[key]);
              }
            }
          }
          decompressed = petri::addExact(decompressed,represents);
        } catch (std::overflow_error &e) {
          std::cerr
              << "Overflow of '128 bit unsigned int' when computing number of decompressed invariant."
              << std::endl;
        }
      }

      return decompressed;
    }



    static void printInvariant (const MatrixCol<T> &invariants,
                                const Permutations &permutations,
                                const std::vector<std::string> &pnames,
                                const std::vector<T> &initial)
    {
      printInvariant (invariants, permutations, pnames, initial, std::cout);
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

    static Invariants computePInvariants (const MatrixCol<T> &pn)
    {
      return computePInvariants (pn, false, 120);
    }

    static Invariants computePInvariants (
        MatrixCol<T> &pn, bool onlyPositive, int timeout,
        const EliminationHeuristic &heuristic = EliminationHeuristic ())
    {
      std::promise<Invariants> promise;
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
        return {MatrixCol<T> (),{}};
      }
    }

  public:

    static Invariants computePInvariants (
        MatrixCol<T> &pn, bool onlyPositive,
        const EliminationHeuristic &heuristic = EliminationHeuristic ())
    {
      auto startTime = std::chrono::steady_clock::now ();
      try {
        auto tpn = pn.transpose ();
        auto inv = petri::InvariantCalculator<T>::calcInvariantsPIPE (
            tpn, onlyPositive, heuristic);
        std::string logMessage = "Computed "
            + std::to_string (inv.first.getColumnCount ()) + " invariants in "
            + std::to_string (
                std::chrono::duration_cast<std::chrono::milliseconds> (
                    std::chrono::steady_clock::now () - startTime).count ())
            + " ms";
        writeToLog (logMessage);
        return inv;
      } catch (std::overflow_error &e) {
        std::cerr << e.what () << std::endl;
        std::string logMessage = "Invariants computation overflowed in "
            + std::to_string (
                std::chrono::duration_cast<std::chrono::milliseconds> (
                    std::chrono::steady_clock::now () - startTime).count ())
            + " ms";
        writeToLog (logMessage);
        return {MatrixCol<T> (),{}};
      }
    }

  };

} /* namespace petri */

#endif /* INVARIANTMIDDLE_H_ */
