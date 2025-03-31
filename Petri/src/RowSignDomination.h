#ifndef ROWSIGNS_DOM_H_
#define ROWSIGNS_DOM_H_

#include <vector>
#include <cstddef>
#include <iosfwd>
#include <map>
#include "SparseArray.h"
#include "MatrixCol.h"
#include "SparseBoolArray.h"


/// Holds sign information for a single row, distinguishing basis (pure positive)
/// and non-basis (mixed-sign or non-basis) columns.
template<typename T>
class RowSignBasis {
public:
  struct PlusSigns {
    SparseBoolArray basis;    // Pure positive basis columns
    SparseBoolArray nonBasis; // Mixed-sign or non-basis columns

    // Combined size for compatibility with existing iteration patterns
    size_t size() const { return basis.size() + nonBasis.size(); }

    // Combined key access for legacy iteration
    size_t keyAt(size_t i) const {
      if (i < basis.size()) return basis.keyAt(i);
      return nonBasis.keyAt(i - basis.size());
    }

    // Access to basis-only columns for minimal support checks
    const SparseBoolArray& basisOnly() const { return basis; }
  };

  const size_t row;         // Row index
  PlusSigns pPlus;          // Positive columns (split into basis and non-basis)
  SparseBoolArray pMinus;   // Negative columns

  /// Constructs a RowSignBasis for the given row, initially empty.
  RowSignBasis(size_t row)
      : row(row), pPlus(), pMinus() {}

  /// Updates the sign information for column 'j' with value 'val'.
  /// If 'inBasis' is true, only modifies pPlus.basis; otherwise, handles
  /// pMinus and pPlus.nonBasis like original RowSign.
  void setValue(size_t j, T val, bool inBasis = false) {
    if (inBasis) {
      if (val == 0) {
        pPlus.basis.clear(j); // Only clear basis if zero
      } else if (val > 0) {
        pPlus.basis.set(j);   // Set in basis, no clears elsewhere
      }
      // No action for val < 0, as basis is pure positive
    } else {
      if (val == 0) {
        pMinus.clear(j);
        pPlus.nonBasis.clear(j);
      } else if (val < 0) {
        pMinus.set(j);
        pPlus.nonBasis.clear(j);
      } else { // val > 0
        pMinus.clear(j);
        pPlus.nonBasis.set(j);
      }
    }
  }

  /// Prints the RowSignBasis to the given stream.
  void print(std::ostream& os) const {
    os << "RowSignBasis [row=" << row << ", pPlus={basis={";
    for (size_t i = 0; i < pPlus.basis.size(); ++i) {
      os << (i > 0 ? ", " : "") << pPlus.basis.keyAt(i);
    }
    os << "}, nonBasis={";
    for (size_t i = 0; i < pPlus.nonBasis.size(); ++i) {
      os << (i > 0 ? ", " : "") << pPlus.nonBasis.keyAt(i);
    }
    os << "}}, pMinus={";
    for (size_t i = 0; i < pMinus.size(); ++i) {
      os << (i > 0 ? ", " : "") << pMinus.keyAt(i);
    }
    os << "}]";
  }

  friend std::ostream& operator<<(std::ostream& os, const RowSignBasis& rs) {
    rs.print(os);
    return os;
  }
};


template<typename T>
class RowSignsDomination {
public:
  using MapType = std::unordered_map<size_t, RowSignBasis<T>>;

private:
  MapType rows;                    // All non-zero rows
  std::unordered_set<size_t> rowsWithSomeNeg; // Rows with any negative entries

  void updateRowsWithSomeNeg(size_t row, const RowSignBasis<T>& rs) {
    bool hasNeg = rs.pMinus.size() > 0;
    if (hasNeg) {
      rowsWithSomeNeg.insert(row);
    } else {
      rowsWithSomeNeg.erase(row);
    }
  }

public:
  class const_value_iterator {
    const std::unordered_set<size_t>* negRows;
    const MapType* rowsMap;
    typename std::unordered_set<size_t>::const_iterator it;
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = const RowSignBasis<T>;
    using difference_type = std::ptrdiff_t;
    using pointer = const RowSignBasis<T>*;
    using reference = const RowSignBasis<T>&;

    const_value_iterator(
        const std::unordered_set<size_t>* nr, const MapType* rm,
        typename std::unordered_set<size_t>::const_iterator i)
        : negRows(nr), rowsMap(rm), it(i) {}
    const_value_iterator& operator++() {
      ++it;
      return *this;
    }
    const_value_iterator operator++(int) {
      auto tmp = *this;
      ++it;
      return tmp;
    }
    reference operator*() const {
      return rowsMap->at(*it);
    }
    pointer operator->() const {
      return &rowsMap->at(*it);
    }
    bool operator==(const const_value_iterator& other) const {
      return it == other.it;
    }
    bool operator!=(const const_value_iterator& other) const {
      return it != other.it;
    }
  };

  /// Constructs from MatrixCol<T>, initializing all purePositive columns as basis.
  RowSignsDomination(const MatrixCol<T>& matC) {
    rows.reserve(matC.getRowCount());
    rowsWithSomeNeg.reserve(matC.getRowCount());

    for (size_t icol = 0, cole = matC.getColumnCount(); icol < cole; ++icol) {
      const SparseArray<T>& col = matC.getColumn(icol);
      bool isBasis = col.isPurePositive(); // Check once per column
      for (size_t i = 0, ie = col.size(); i < ie; ++i) {
        size_t row = col.keyAt(i);
        T val = col.valueAt(i);
        auto [it, inserted] = rows.emplace(row, RowSignBasis<T>(row));
        RowSignBasis<T>& rs = it->second;
        if (isBasis) {
          rs.pPlus.basis.append(icol, true); // Pure positive, always positive
        } else if (val < 0) {
          rs.pMinus.append(icol, true);
        } else if (val > 0) {
          rs.pPlus.nonBasis.append(icol, true);
        }
        // No action for val == 0, as we’re building from non-zero entries
      }
    }

    // Populate rowsWithSomeNeg after building rows
    for (const auto& [row, rs] : rows) {
      if (rs.pMinus.size() > 0) {
        rowsWithSomeNeg.insert(row);
      }
    }
  }

  void setValue(size_t row, size_t col, T newVal, bool inBasis = false) {
    auto [it, inserted] = rows.emplace(row, RowSignBasis<T>(row));
    RowSignBasis<T>& rs = it->second;
    rs.setValue(col, newVal, inBasis);

    if (rs.pPlus.size() == 0 && rs.pMinus.size() == 0) {
      rows.erase(it);
      rowsWithSomeNeg.erase(row);
    } else {
      updateRowsWithSomeNeg(row, rs);
    }
  }

  const RowSignBasis<T>& get(size_t row) const {
    auto it = rows.find(row);
    if (it != rows.end()) {
      return it->second;
    }
    static const RowSignBasis<T> empty(row);
    return empty;
  }

  size_t size() const {
    return rows.size();
  }

  const_value_iterator begin() const {
    return const_value_iterator(&rowsWithSomeNeg, &rows, rowsWithSomeNeg.begin());
  }
  const_value_iterator end() const {
    return const_value_iterator(&rowsWithSomeNeg, &rows, rowsWithSomeNeg.end());
  }

  ssize_t findSingleSignRow(size_t limit) const {
    size_t count = 0;
    for (const auto& rowIdx : rowsWithSomeNeg) {
      const RowSignBasis<T>& rs = rows.at(rowIdx);
      if (rs.pPlus.size() == 1 || rs.pMinus.size() == 1) {
        return static_cast<ssize_t>(rowIdx);
      }
      if (++count >= limit) break;
    }
    return -1;
  }

  void clearRow(size_t row) {
    rows.erase(row);
    rowsWithSomeNeg.erase(row);
  }

  friend std::ostream& operator<<(std::ostream& os, const RowSignsDomination<T>& srs) {
    os << "RowSignsBasisDomination (" << srs.rows.size() << " total rows, "
       << srs.rowsWithSomeNeg.size() << " rows with negatives):\n";
    os << "  RowsWithSomeNeg: ";
    for (const auto& idx : srs.rowsWithSomeNeg)
      os << idx << " ";
    os << "\n  All Rows:\n";
    for (const auto& [idx, rs] : srs.rows)
      os << rs << "\n";
    return os;
  }
};


#endif /* ROWSIGNS_DOM_H_ */
