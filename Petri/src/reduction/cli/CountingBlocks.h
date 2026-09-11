#pragma once
#include <ostream>
#include <string>
#include "io/PNETIO.h"
#include "reduction/Counting.h"

namespace petri::reduction {
/** The counting record as PNET named blocks and back (io/PNET.md). A block
 * is one column; `TMULT` and `PCOEF` are indexed by object, `PDROP` is a list.
 * Reading takes the three names and reports the others, which nothing here
 * maintains; a `TMULT` is required for the arcs to count as vouched for. */
template<class T>
Counting<T> countingFromBlocks(const typename PNETIO<T>::Blocks& blocks, size_t places,
                               size_t transitions, std::ostream& diagnostics) {
  Counting<T> c;
  c.pcoef.assign(places, T(0));
  c.arcsLost = "the input net carries no TMULT block";
  for (const auto& [name, matrix] : blocks) {
    if (matrix.getColumnCount() != 1) throw std::invalid_argument("PNET block " + name + " is not one column");
    const auto& col = matrix.getColumn(0);
    if (name == "TMULT") {
      if (matrix.getRowCount() != transitions) throw std::invalid_argument("TMULT rows differ from the transition count");
      c.tmult = std::vector<T>(transitions, T(0)); c.arcsLost.clear();
      for (size_t i = 0; i < col.size(); ++i) (*c.tmult)[col.keyAt(i)] = col.valueAt(i);
    } else if (name == "PCOEF") {
      if (matrix.getRowCount() != places) throw std::invalid_argument("PCOEF rows differ from the place count");
      for (size_t i = 0; i < col.size(); ++i) c.pcoef[col.keyAt(i)] = col.valueAt(i);
    } else if (name == "PDROP") {
      for (size_t i = 0; i < col.size(); ++i) c.pdrop.push_back(col.valueAt(i));
    } else {
      diagnostics << "Reduction drops the input block " << name << ": nothing maintains it.\n";
    }
  }
  return c;
}

/** `TMULT` while the arcs are vouched for, `PDROP` when constant places
 * holding tokens went, `PCOEF` when some place stands for several. */
template<class T>
typename PNETIO<T>::Blocks countingToBlocks(const Counting<T>& c) {
  typename PNETIO<T>::Blocks blocks;
  auto column = [](const std::vector<T>& values, bool asList) {
    MatrixCol<T> m(values.size(), 0);
    SparseArray<T> col;
    for (size_t i = 0; i < values.size(); ++i) if (asList || values[i] != 0) col.append(i, values[i]);
    m.appendColumn(std::move(col));
    return m;
  };
  if (c.tmult) blocks.emplace_back("TMULT", column(*c.tmult, false));
  if (!c.pdrop.empty()) blocks.emplace_back("PDROP", column(c.pdrop, true));
  if (c.weighted()) blocks.emplace_back("PCOEF", column(c.pcoef, false));
  return blocks;
}

/** One line per fact of the record, for the diagnostics. */
template<class T>
void describeCounting(const Counting<T>& c, std::ostream& os) {
  if (c.tmult) {
    size_t fused = 0; T extra = 0;
    for (T m : *c.tmult) if (m != 0) { ++fused; extra = petri::addExact(extra, m); }
    os << "Counting record: TMULT kept, " << fused << " transitions stand for " << extra << " more.\n";
  } else os << "Counting record: TMULT dropped (" << c.arcsLost << ").\n";
  if (!c.pdrop.empty()) {
    T total = 0; for (T m : c.pdrop) total = petri::addExact(total, m);
    os << "Counting record: PDROP " << c.pdrop.size() << " constant places held " << total << " tokens.\n";
  }
  if (c.weighted()) {
    size_t fused = 0; T extra = 0;
    for (T k : c.pcoef) if (k != 0) { ++fused; extra = petri::addExact(extra, k); }
    os << "Counting record: PCOEF " << fused << " places stand for " << extra << " more.\n";
  }
}
}
