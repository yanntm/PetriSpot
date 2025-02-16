#ifndef ROWSIGNS_H_
#define ROWSIGNS_H_

#include <vector>
#include <cstddef>
#include <iosfwd>
#include "SparseArray.h"
#include "MatrixCol.h"
#include "SparseBoolArray.h"

namespace petri {

/// Holds the sign information for a single row in the incidence matrix.
template<typename T>
class RowSign {
public:
  // The row index.
  const size_t row;
  // Set of column indices where the entry is positive.
  SparseBoolArray pPlus;
  // Set of column indices where the entry is negative.
  SparseBoolArray pMinus;

  /// Constructs a RowSign for the given row, initially empty.
  RowSign(size_t row)
    : row(row), pPlus(), pMinus()
  { }

  /// Updates the sign information for column 'j' with value 'val'.
  void setValue(size_t j, T val) {
    if (val == 0) {
      pMinus.clear(j);
      pPlus.clear(j);
    } else if (val < 0) {
      pMinus.set(j);
      pPlus.clear(j);
    } else {
      pMinus.clear(j);
      pPlus.set(j);
    }
  }

  /// Prints the RowSign to the given stream.
  void print(std::ostream &os) const {
    os << "RowSign [row=" << row << ", pPlus=" << pPlus << ", pMinus=" << pMinus << "]";
  }

  friend std::ostream& operator<<(std::ostream &os, const RowSign &rs) {
    rs.print(os);
    return os;
  }
};

/// Holds the sign information for each row in the incidence matrix.
template<typename T>
class RowSigns {
  // The container holding a RowSign for each row.
  std::vector<RowSign<T>> rows;

public:
  /// Constructs the RowSigns data from the incidence matrix 'matC'.
  /// For each nonzero entry in 'matC', the corresponding row's sign set is updated.
  RowSigns(const MatrixCol<T> &matC) {
    size_t rowCount = matC.getRowCount();
    rows.reserve(rowCount);
    for (size_t row = 0; row < rowCount; ++row) {
      rows.emplace_back(row);
    }
    for (size_t icol = 0, cole = matC.getColumnCount(); icol < cole; ++icol) {
      const SparseArray<T> &col = matC.getColumn(icol);
      for (size_t i = 0, ie = col.size(); i < ie; ++i) {
        size_t row = col.keyAt(i);
        if (col.valueAt(i) < 0) {
          rows[row].pMinus.append(icol, true);
        } else {
          rows[row].pPlus.append(icol, true);
        }
      }
    }
  }

  /// Forces all modifications to row sign data to go through this method.
  void setValue(size_t row, size_t col, T newVal) {
    rows[row].setValue(col, newVal);
  }

  /// Returns a const reference to the RowSign for the given row.
  const RowSign<T>& get(size_t row) const {
    return rows[row];
  }

  /// Returns the number of rows.
  size_t size() const {
    return rows.size();
  }

  /// Scans from startIndex (wrapping around) to find the first row where either
  /// pPlus or pMinus has exactly one element. Returns the row index or -1 if none found.
  ssize_t findSingleSignRow(size_t startIndex) const {
    size_t sz = rows.size();
    for (size_t i = startIndex; i < sz; ++i) {
      if (rows[i].pMinus.size() == 1 || rows[i].pPlus.size() == 1)
        return static_cast<ssize_t>(i);
    }
    for (size_t i = 0; i < startIndex; ++i) {
      if (rows[i].pMinus.size() == 1 || rows[i].pPlus.size() == 1)
        return static_cast<ssize_t>(i);
    }
    return -1;
  }
};

} // namespace petri

#endif // ROWSIGNS_H_
