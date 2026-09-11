#pragma once

#include <chrono>
#include <limits>
#include <optional>
#include <string>
#include <vector>
#include "core/SparsePetriNet.h"
#include "reduction/Configuration.h"

namespace petri::reduction {

/** Sparse edit workspace. Slots stay stable until publication. Every adjacency
 * row contains active arcs only; inactive empty columns are never transitions.
 * Observed places are retained, making output remapping a pure permutation. */
template<class T> class Workspace {
  void replace(MatrixCol<T>& matrix, MatrixCol<T>& transpose, size_t t,
               SparseArray<T> column) {
    const auto& old = matrix.getColumn(t);
    for (size_t i = 0; i < old.size(); ++i)
      transpose.getColumn(old.keyAt(i)).put(t, 0);
    for (size_t i = 0; i < column.size(); ++i)
      transpose.getColumn(column.keyAt(i)).put(t, column.valueAt(i));
    matrix.getColumn(t) = std::move(column);
  }

public:
  MatrixCol<T> pre, post, consumers, producers;
  std::vector<T> marks;
  std::vector<std::string> places, transitions;
  std::string name;
  std::vector<bool> liveP, liveT, observed;
  Configuration config;
  size_t changes = 0;
  bool limited = false;
  bool safe = false; // the input's one-safety; a rule that fuses places clears it
  std::optional<bool> deadlock;
  std::chrono::steady_clock::time_point deadline;

  Workspace(SparsePetriNet<T> net, Configuration options, std::vector<bool> support)
      : pre(std::move(net.getFlowPT())), post(std::move(net.getFlowTP())),
        marks(net.getMarks()), places(net.getPnames()), transitions(net.getTnames()),
        name(net.getName()), liveP(places.size(), true), liveT(transitions.size(), true),
        observed(std::move(support)), config(options), safe(net.isSafe()),
        deadline(std::chrono::steady_clock::now() + options.timeLimit) {
    if (observed.empty()) observed.resize(places.size(), false);
    if (observed.size() != places.size() || pre.getRowCount() != places.size()
        || post.getRowCount() != places.size() || pre.getColumnCount() != transitions.size()
        || post.getColumnCount() != transitions.size())
      throw std::invalid_argument("Inconsistent reduction net dimensions/support");
    for (T m : marks) if (m < 0) throw std::invalid_argument("Negative marking");
    for (const auto* mat : {&pre, &post})
      for (size_t t = 0; t < transitions.size(); ++t) {
        const auto& col = mat->getColumn(t);
        for (size_t i = 0; i < col.size(); ++i)
          if (col.keyAt(i) >= places.size() || col.valueAt(i) <= 0)
            throw std::invalid_argument("Invalid reduction arc");
      }
    consumers = pre.transpose();
    producers = post.transpose();
  }

  bool stop() {
    if (config.timeLimit.count() > 0 && std::chrono::steady_clock::now() >= deadline)
      limited = true;
    return limited;
  }

  bool visible(size_t t) const {
    for (const auto* mat : {&pre, &post}) {
      const auto& col = mat->getColumn(t);
      for (size_t i = 0; i < col.size(); ++i) {
        size_t p = col.keyAt(i);
        if (observed[p] && pre.getColumn(t).get(p) != post.getColumn(t).get(p)) return true;
      }
    }
    return false;
  }

  void replacePre(size_t t, SparseArray<T> col) {
    replace(pre, consumers, t, std::move(col)); ++changes;
  }
  void replacePost(size_t t, SparseArray<T> col) {
    replace(post, producers, t, std::move(col)); ++changes;
  }
  size_t appendTransition(SparseArray<T> input, SparseArray<T> output, std::string label) {
    size_t t = transitions.size();
    pre.appendColumn(SparseArray<T>{}); post.appendColumn(SparseArray<T>{});
    consumers.addRow(); producers.addRow();
    transitions.push_back(std::move(label)); liveT.push_back(true);
    replacePre(t, std::move(input)); replacePost(t, std::move(output));
    return t;
  }

  void retireTransition(size_t t) {
    if (!liveT[t]) return;
    replacePre(t, {}); replacePost(t, {}); liveT[t] = false;
  }

  /** Remove a proven irrelevant guard/coordinate; no transition is implicitly
   * retired. Callers first retire consumers proved dead by this coordinate. */
  void erasePlaceArcs(size_t p) {
    for (size_t i = 0; i < consumers.getColumn(p).size(); ++i)
      pre.getColumn(consumers.getColumn(p).keyAt(i)).put(p, 0);
    for (size_t i = 0; i < producers.getColumn(p).size(); ++i)
      post.getColumn(producers.getColumn(p).keyAt(i)).put(p, 0);
    if (consumers.getColumn(p).size() || producers.getColumn(p).size()) ++changes;
    consumers.getColumn(p).clear(); producers.getColumn(p).clear();
  }
  void retirePlace(size_t p) {
    if (!liveP[p] || observed[p]) throw std::logic_error("Retiring a protected/inactive place");
    erasePlaceArcs(p); liveP[p] = false; ++changes;
  }

  std::string composedName(size_t h, size_t f) const {
    return transitions[h] + "." + transitions[f];
  }

  SparsePetriNet<T> publish(std::vector<size_t>& placeMap,
                           std::vector<size_t>& transitionMap) const {
    constexpr size_t absent = std::numeric_limits<size_t>::max();
    placeMap.assign(places.size(), absent); transitionMap.assign(transitions.size(), absent);
    SparsePetriNet<T> result;
    result.setName(name);
    result.setSafe(safe);
    for (size_t p = 0; p < places.size(); ++p)
      if (liveP[p]) placeMap[p] = result.addPlace(places[p], marks[p]);
    for (size_t t = 0; t < transitions.size(); ++t) if (liveT[t]) {
      size_t out = result.addTransition(transitions[t]); transitionMap[t] = out;
      for (const auto* mat : {&pre, &post}) {
        const auto& col = mat->getColumn(t);
        for (size_t i = 0; i < col.size(); ++i) {
          size_t p = placeMap[col.keyAt(i)];
          if (p == absent) throw std::logic_error("Arc to retired place");
          if (p > static_cast<size_t>(std::numeric_limits<int>::max())
              || out > static_cast<size_t>(std::numeric_limits<int>::max()))
            throw std::overflow_error("Net builder index exceeds int");
          if (mat == &pre) result.addPreArc(static_cast<int>(p), static_cast<int>(out), col.valueAt(i));
          else result.addPostArc(static_cast<int>(p), static_cast<int>(out), col.valueAt(i));
        }
      }
    }
    return result;
  }
};
} // namespace petri::reduction
