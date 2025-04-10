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

  static void confrontIndex(MatrixCol<T>& colsB, RowSignsDomination<T>& rowSigns, std::unordered_set<size_t>& basisIndices) {
    // Step 1: Verify basisIndices against pure positive columns in colsB
    std::unordered_set<size_t> trueBasisIndices;
    for (size_t icol = 0; icol < colsB.getColumnCount(); ++icol) {
      const SparseArray<T>& col = colsB.getColumn(icol);
      if (col.isPurePositive() && col.size() > 0) {
        trueBasisIndices.insert(icol);
      }
    }

    // Check for columns in basisIndices that are not pure positive
    bool basisMismatch = false;
    for (size_t icol : basisIndices) {
      if (trueBasisIndices.find(icol) == trueBasisIndices.end()) {
        std::cout << "Error: Column " << icol << " is in basisIndices but not pure positive.\n";
        std::cout << "Column " << icol << " contents: " << colsB.getColumn(icol) << "\n";
        basisMismatch = true;
      }
    }

    // Check for pure positive columns not in basisIndices
    for (size_t icol : trueBasisIndices) {
      if (basisIndices.find(icol) == basisIndices.end()) {
        std::cout << "Error: Column " << icol << " is pure positive but not in basisIndices.\n";
        std::cout << "Column " << icol << " contents: " << colsB.getColumn(icol) << "\n";
        basisMismatch = true;
      }
    }

    if (!basisMismatch) {
      std::cout << "basisIndices is consistent with pure positive columns in colsB.\n";
    } else {
      std::cout << "colsB: " << colsB << "\n";
      std::cout << "basisIndices: { ";
      for (size_t idx : basisIndices) std::cout << idx << " ";
      std::cout << "}\n";
    }

    // Step 2: Construct a fresh RowSignsDomination<T> using colsB and basisIndices
    RowSignsDomination<T> freshRowSigns(colsB);

    // Step 3: Compare freshRowSigns with existing rowSigns
    bool rowSignsMismatch = false;

    // Get row sets
    std::unordered_set<size_t> freshRows, existingRows;
    for (const auto& [row, _] : freshRowSigns.rows) freshRows.insert(row);
    for (const auto& [row, _] : rowSigns.rows) existingRows.insert(row);

    if (freshRows != existingRows) {
      std::cout << "Mismatch in row sets between rowSigns and freshRowSigns:\n";
      for (size_t row : freshRows) {
        if (existingRows.count(row) == 0) {
          std::cout << "Row " << row << " in freshRowSigns but not in rowSigns.\n";
        }
      }
      for (size_t row : existingRows) {
        if (freshRows.count(row) == 0) {
          std::cout << "Row " << row << " in rowSigns but not in freshRowSigns.\n";
        }
      }
      rowSignsMismatch = true;
    }

    // Compare each row's sign information
    for (size_t row : freshRows) {
      if (existingRows.count(row) == 0) continue; // Skip rows not in both
      const RowSignBasis<T>& freshRS = freshRowSigns.get(row);
      const RowSignBasis<T>& existingRS = rowSigns.get(row);

      if (freshRS.pPlus.basis != existingRS.pPlus.basis) {
        std::cout << "Mismatch in pPlus.basis for row " << row << ":\n";
        std::cout << "Fresh: " << freshRS.pPlus.basis << "\n";
        std::cout << "Existing: " << existingRS.pPlus.basis << "\n";
        rowSignsMismatch = true;
      }
      if (freshRS.pPlus.nonBasis != existingRS.pPlus.nonBasis) {
        std::cout << "Mismatch in pPlus.nonBasis for row " << row << ":\n";
        std::cout << "Fresh: " << freshRS.pPlus.nonBasis << "\n";
        std::cout << "Existing: " << existingRS.pPlus.nonBasis << "\n";
        rowSignsMismatch = true;
      }
      if (freshRS.pMinus != existingRS.pMinus) {
        std::cout << "Mismatch in pMinus for row " << row << ":\n";
        std::cout << "Fresh: " << freshRS.pMinus << "\n";
        std::cout << "Existing: " << existingRS.pMinus << "\n";
        rowSignsMismatch = true;
      }
    }

    // Compare rowsWithSomeNeg
    std::unordered_set<size_t> freshRowsWithSomeNeg;
    for (const auto& [row, rs] : freshRowSigns.rows) {
      if (rs.pMinus.size() > 0) freshRowsWithSomeNeg.insert(row);
    }
    if (freshRowsWithSomeNeg != rowSigns.rowsWithSomeNeg) {
      std::cout << "Mismatch in rowsWithSomeNeg:\n";
      for (size_t row : freshRowsWithSomeNeg) {
        if (rowSigns.rowsWithSomeNeg.count(row) == 0) {
          std::cout << "Row " << row << " has negatives in freshRowSigns but not in rowSigns.\n";
        }
      }
      for (size_t row : rowSigns.rowsWithSomeNeg) {
        if (freshRowsWithSomeNeg.count(row) == 0) {
          std::cout << "Row " << row << " has negatives in rowSigns but not in freshRowSigns.\n";
        }
      }
      rowSignsMismatch = true;
    }

    if (!rowSignsMismatch) {
      std::cout << "rowSigns is consistent with colsB and basisIndices.\n";
    } else {
      std::cout << "Full rowSigns: " << rowSigns << "\n";
      std::cout << "Full freshRowSigns: " << freshRowSigns << "\n";
    }
  }
};


#endif /* ROWSIGNS_DOM_H_ */
