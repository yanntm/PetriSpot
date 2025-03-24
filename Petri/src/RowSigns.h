#ifndef ROWSIGNS_H_
#define ROWSIGNS_H_

#include <vector>
#include <cstddef>
#include <iosfwd>
#include <map>
#include "SparseArray.h"
#include "MatrixCol.h"
#include "SparseBoolArray.h"

namespace petri
{

/// Holds the sign information for a single row in the incidence matrix.
template<typename T>
  class RowSign
  {
  public:
    // The row index.
    const size_t row;
    // Set of column indices where the entry is positive.
    SparseBoolArray pPlus;
    // Set of column indices where the entry is negative.
    SparseBoolArray pMinus;

    /// Constructs a RowSign for the given row, initially empty.
    RowSign (size_t row)
        : row (row), pPlus (), pMinus ()
    {
    }

    /// Updates the sign information for column 'j' with value 'val'.
    void setValue (size_t j, T val)
    {
      if (val == 0) {
        pMinus.clear (j);
        pPlus.clear (j);
      } else if (val < 0) {
        pMinus.set (j);
        pPlus.clear (j);
      } else {
        pMinus.clear (j);
        pPlus.set (j);
      }
    }

    /// Prints the RowSign to the given stream.
    void print (std::ostream &os) const
    {
      os << "RowSign [row=" << row << ", pPlus=" << pPlus << ", pMinus="
          << pMinus << "]";
    }

    friend std::ostream& operator<< (std::ostream &os, const RowSign &rs)
    {
      rs.print (os);
      return os;
    }
  };

/// Holds the sign information for each row in the incidence matrix.
template<typename T>
  class DenseRowSigns
  {
    // The container holding a RowSign for each row.
    using RowsType = std::vector<RowSign<T>>;
    RowsType rows;

  public:
    /// Constructs the RowSigns data from the incidence matrix 'matC'.
    /// For each nonzero entry in 'matC', the corresponding row's sign set is updated.
    DenseRowSigns (const MatrixCol<T> &matC)
    {
      size_t rowCount = matC.getRowCount ();
      rows.reserve (rowCount);
      for (size_t row = 0; row < rowCount; ++row) {
        rows.emplace_back (row);
      }
      for (size_t icol = 0, cole = matC.getColumnCount (); icol < cole;
          ++icol) {
        const SparseArray<T> &col = matC.getColumn (icol);
        for (size_t i = 0, ie = col.size (); i < ie; ++i) {
          size_t row = col.keyAt (i);
          if (col.valueAt (i) < 0) {
            rows[row].pMinus.append (icol, true);
          } else {
            rows[row].pPlus.append (icol, true);
          }
        }
      }
    }

    /// Forces all modifications to row sign data to go through this method.
    void setValue (size_t row, size_t col, T newVal)
    {
      rows[row].setValue (col, newVal);
    }

    /// Returns a const reference to the RowSign for the given row.
    const RowSign<T>& get (size_t row) const
    {
      return rows[row];
    }

    /// Returns the number of rows.
    size_t size () const
    {
      return rows.size ();
    }

    /// Returns a const iterator to the first stored RowSign.
    typename RowsType::const_iterator begin () const
    {
      return rows.begin ();
    }

    /// Returns a const iterator to one past the last stored RowSign.
    typename RowsType::const_iterator end () const
    {
      return rows.end ();
    }

    /// Scans from startIndex (wrapping around) to find the first row where either
    /// pPlus or pMinus has exactly one element. Returns the row index or -1 if none found.
    ssize_t findSingleSignRow (size_t startIndex) const
    {
      size_t sz = rows.size ();
      for (size_t i = startIndex; i < sz; ++i) {
        if (rows[i].pMinus.size () == 1 || rows[i].pPlus.size () == 1) return static_cast<ssize_t> (i);
      }
      for (size_t i = 0; i < startIndex; ++i) {
        if (rows[i].pMinus.size () == 1 || rows[i].pPlus.size () == 1) return static_cast<ssize_t> (i);
      }
      return -1;
    }

//  ssize_t findSingleSignRow () const
//    {
//      // Iterate only over stored (nonempty) row signs.
//      for (const auto &rs : *this) { // assuming RowSigns supports begin() and end()
//        if (rs.pPlus.size () == 1 || rs.pMinus.size () == 1) return static_cast<ssize_t> (rs.row);
//      }
//      return -1;
//    }
  };

template<typename T>
  class SparseRowSigns
  {
  public:
    using MapType = std::unordered_map<size_t, RowSign<T>>; // Key is row index.

  private:
    MapType rows;
    bool indexSingleSignRows;
    bool ignorePurePositiveRows;
    mutable std::vector<size_t> singleSignRows;
    std::unordered_set<size_t> purePositiveRows;
  public:
    // --- Iterator API ---
    // A const iterator that yields a const reference to a RowSign.
    class const_value_iterator
    {
      typename MapType::const_iterator it;
    public:
      using iterator_category = std::forward_iterator_tag;
      using value_type = const RowSign<T>;
      using difference_type = std::ptrdiff_t;
      using pointer = const RowSign<T>*;
      using reference = const RowSign<T>&;

      const_value_iterator (typename MapType::const_iterator it)
          : it (it)
      {
      }
      const_value_iterator& operator++ ()
      {
        ++it;
        return *this;
      }
      const_value_iterator operator++ (int)
      {
        const_value_iterator tmp (*this);
        ++it;
        return tmp;
      }
      reference operator* () const
      {
        return it->second;
      }
      pointer operator-> () const
      {
        return &(it->second);
      }
      bool operator== (const const_value_iterator &other) const
      {
        return it == other.it;
      }
      bool operator!= (const const_value_iterator &other) const
      {
        return it != other.it;
      }
    };

  public:
    // --- Constructors ---
    /// Constructs the SparseRowSigns from the incidence matrix 'matC'
    /// using a two‑pass approach: first build a dense RowSigns,
    /// then insert only those rows with nonzero sign data into the map.
    SparseRowSigns (const MatrixCol<T> &matC, bool indexSingleSignRows = false,
                    bool ignorePurePos = false)
        : indexSingleSignRows (indexSingleSignRows), ignorePurePositiveRows (
            ignorePurePos)
    {
      // Build the dense version.
      DenseRowSigns<T> dense (matC);
      rows.reserve (matC.getRowCount ());
      // Insert only nonzero rows.
      for (size_t i = 0; i < dense.size (); i++) {
        const auto &rs = dense.get (i);
        if (ignorePurePos && rs.pPlus.size () != 0 && rs.pMinus.size () == 0) {
          purePositiveRows.insert (rs.row);
          continue;
        }
        if (rs.pPlus.size () != 0 || rs.pMinus.size () != 0) {
          rows.insert (
            { rs.row, rs });
        }
        if (indexSingleSignRows
            && (rs.pPlus.size () == 1 || rs.pMinus.size () == 1)) {
          singleSignRows.push_back (rs.row);
        }
      }
    }

    // --- Modification ---
    /// Forces all modifications to row sign data to go through this method.
    /// If the update causes the RowSign to become empty, erase the entry.
    void setValue (size_t row, size_t col, T newVal)
    {
      if (ignorePurePositiveRows
          && purePositiveRows.find (row) != purePositiveRows.end ()) return;
      auto it = rows.find (row);
      if (it == rows.end ()) {
        if (newVal != 0) {
          if (ignorePurePositiveRows && newVal > 0) {
            purePositiveRows.insert (row);
            return;
          }
          RowSign<T> rs (row);
          rs.setValue (col, newVal);
          rows.insert (
            { row, std::move (rs) });
          if (indexSingleSignRows) {
            singleSignRows.push_back (rs.row);
          }
        }
      } else {
        it->second.setValue (col, newVal);
        // Only remove if newVal is zero and both sets become empty.
        if (newVal == 0 && it->second.pPlus.size () == 0
            && it->second.pMinus.size () == 0) {
          rows.erase (it);
        } else if (ignorePurePositiveRows && newVal >= 0
            && it->second.pMinus.size () == 0) {
          purePositiveRows.insert (row);
          rows.erase (it);
          return;
        } else if (indexSingleSignRows
            && (it->second.pPlus.size () == 1 || it->second.pMinus.size () == 1)) {
          singleSignRows.push_back (it->second.row);
        }
      }
    }

    // --- Accessors ---
    /// Returns a const reference to the RowSign for the given row.
    /// If the row is not stored, returns a static default empty RowSign.
    const RowSign<T>& get (size_t row) const
    {
      auto it = rows.find (row);
      if (it != rows.end ()) return it->second;
      static const RowSign<T> empty (row);
      return empty;
    }

    /// Returns the number of stored (nonzero) RowSign entries.
    size_t size () const
    {
      return rows.size ();
    }

    /// Returns a const iterator to the first stored RowSign.
    const_value_iterator begin () const
    {
      return const_value_iterator (rows.begin ());
    }

    /// Returns a const iterator to one past the last stored RowSign.
    const_value_iterator end () const
    {
      return const_value_iterator (rows.end ());
    }

    // --- Candidate Search ---
    /// Iterates over the stored rows (in order) to find the first row where either
    /// pPlus or pMinus has exactly one element. Returns that row index, or -1 if none found.
    ssize_t findSingleSignRow (size_t limit) const
    {
      if (indexSingleSignRows) {
        for (size_t index = singleSignRows.size (); index-- > 0;) {
          const auto &rs = get (singleSignRows[index]);
          singleSignRows.pop_back ();
          if (rs.pPlus.size () == 1 || rs.pMinus.size () == 1) return static_cast<ssize_t> (rs.row);
        }
      } else {
        size_t current = 0;
        for (const auto &rs : *this) {
          if (rs.pPlus.size () == 1 || rs.pMinus.size () == 1) return static_cast<ssize_t> (rs.row);
          if (++current >= limit) break;
        }
      }
      return -1;
    }

    void clearRow (size_t row)
    {
      rows.erase (row);
    }

    friend std::ostream& operator<< (std::ostream &os,
                                     const SparseRowSigns<T> &srs)
    {
      os << "SparseRowSigns (" << srs.rows.size () << " entries):\n";
      os << "  indexSingleSignRows: "
          << (srs.indexSingleSignRows ? "true" : "false") << "\n";
      os << "  ignorePurePositiveRows: "
          << (srs.ignorePurePositiveRows ? "true" : "false") << "\n";
      os << "  Cached singleSignRows: ";
      for (size_t i = 0; i < srs.singleSignRows.size (); ++i)
        os << srs.singleSignRows[i] << " ";
      os << "\n";
      os << "  Cached purePositiveRows: ";
      for (const auto &idx : srs.purePositiveRows)
        os << idx << " ";
      os << "\n";
      os << "  Stored Rows:\n";
      for (auto it = srs.begin (); it != srs.end (); ++it) {
        os << *it << "\n";
      }
      return os;
    }
  };

template<typename T>
  class RowSignsDomination
  {
  public:
    using MapType = std::unordered_map<size_t, RowSign<T>>;

  private:
    MapType rows;                    // All non-zero rows
    std::unordered_set<size_t> rowsWithSomeNeg; // Rows with any negative entries

    void updateRowsWithSomeNeg (size_t row, const RowSign<T> &rs)
    {
      bool hasNeg = rs.pMinus.size () > 0;
      if (hasNeg) {
        rowsWithSomeNeg.insert (row);
      } else {
        rowsWithSomeNeg.erase (row);
      }
    }

  public:
    class const_value_iterator
    {
      const std::unordered_set<size_t> *negRows;
      const MapType *rowsMap;
      typename std::unordered_set<size_t>::const_iterator it;
    public:
      using iterator_category = std::forward_iterator_tag;
      using value_type = const RowSign<T>;
      using difference_type = std::ptrdiff_t;
      using pointer = const RowSign<T>*;
      using reference = const RowSign<T>&;

      const_value_iterator (
          const std::unordered_set<size_t> *nr, const MapType *rm,
          typename std::unordered_set<size_t>::const_iterator i)
          : negRows (nr), rowsMap (rm), it (i)
      {
      }
      const_value_iterator& operator++ ()
      {
        ++it;
        return *this;
      }
      const_value_iterator operator++ (int)
      {
        auto tmp = *this;
        ++it;
        return tmp;
      }
      reference operator* () const
      {
        return rowsMap->at (*it);
      }
      pointer operator-> () const
      {
        return &rowsMap->at (*it);
      }
      bool operator== (const const_value_iterator &other) const
      {
        return it == other.it;
      }
      bool operator!= (const const_value_iterator &other) const
      {
        return it != other.it;
      }
    };

    RowSignsDomination (const MatrixCol<T> &matC)
    {
      // Build dense version
      DenseRowSigns<T> dense (matC);
      rows.reserve (matC.getRowCount ());
      rowsWithSomeNeg.reserve (matC.getRowCount ());

      // Populate sparse rows and rowsWithSomeNeg
      for (size_t i = 0; i < dense.size (); ++i) {
        const auto &rs = dense.get (i);
        if (rs.pPlus.size () != 0 || rs.pMinus.size () != 0) {
          rows.insert (
            { rs.row, rs });
          if (rs.pMinus.size () > 0) {
            rowsWithSomeNeg.insert (rs.row);
          }
        }
      }
    }

    void setValue (size_t row, size_t col, T newVal)
    {
      auto [it, inserted] = rows.emplace (row, RowSign<T> (row));
      RowSign<T> &rs = it->second;
      rs.setValue (col, newVal);

      if (rs.pPlus.size () == 0 && rs.pMinus.size () == 0) {
        rows.erase (it);
        rowsWithSomeNeg.erase (row);
      } else {
        updateRowsWithSomeNeg (row, rs);
      }
    }

    const RowSign<T>& get (size_t row) const
    {
      auto it = rows.find (row);
      if (it != rows.end ()) {
        return it->second;
      }
      static const RowSign<T> empty (row);
      return empty;
    }

    size_t size () const
    {
      return rows.size ();
    }

    const_value_iterator begin () const
    {
      return const_value_iterator (&rowsWithSomeNeg, &rows,
                                   rowsWithSomeNeg.begin ());
    }
    const_value_iterator end () const
    {
      return const_value_iterator (&rowsWithSomeNeg, &rows,
                                   rowsWithSomeNeg.end ());
    }

    ssize_t findSingleSignRow (size_t limit) const
    {
      size_t count = 0;
      for (const auto &rowIdx : rowsWithSomeNeg) {
        const RowSign<T> &rs = rows.at (rowIdx);
        if (rs.pPlus.size () == 1 || rs.pMinus.size () == 1) {
          return static_cast<ssize_t> (rowIdx);
        }
        if (++count >= limit) break;
      }
      return -1;
    }

    void clearRow (size_t row)
    {
      rows.erase (row);
      rowsWithSomeNeg.erase (row);
    }

    friend std::ostream& operator<< (std::ostream &os,
                                     const RowSignsDomination<T> &srs)
    {
      os << "RowSignsDomination (" << srs.rows.size () << " total rows, "
          << srs.rowsWithSomeNeg.size () << " rows with negatives):\n";
      os << "  RowsWithSomeNeg: ";
      for (const auto &idx : srs.rowsWithSomeNeg)
        os << idx << " ";
      os << "\n  All Rows:\n";
      for (const auto& [idx, rs] : srs.rows)
        os << rs << "\n";
      return os;
    }
  };

#define RowSigns SparseRowSigns

} // namespace petri

#endif // ROWSIGNS_H_
