#pragma once

#include <algorithm>
#include <chrono>
#include <limits>
#include <optional>
#include <string>
#include <vector>
#include "core/SparsePetriNet.h"
#include "reduction/Configuration.h"
#include "reduction/Counting.h"

namespace petri::reduction {

/** What the state-equation dead transition tests did, over every pass. */
struct DeadStats {
  size_t found = 0, tested = 0, solves = 0, pivots = 0, passes = 0, cursor = 0;
  size_t places = 0, stuck = 0, byBound = 0; // phase one: places bounded, never marked, transitions dead by a bound
  long placesMs = 0;
  long spentMs = 0;
  bool limited = false;
};

/** Sparse edit workspace. Slots stay stable until publication. Every adjacency
 * row contains active arcs only; inactive empty columns are never transitions.
 * Observed places are retained, making output remapping a pure permutation.
 * `counting`, when attached, is maintained by the operations that say what
 * they did (fuse, dead, constant) and invalidated or refused by the generic
 * ones (Counting.h). */
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
  std::optional<Counting<T>> counting;
  DeadStats dead;
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
    if (counting) throw std::logic_error("Composing a transition under a counting record");
    size_t t = transitions.size();
    pre.appendColumn(SparseArray<T>{}); post.appendColumn(SparseArray<T>{});
    consumers.addRow(); producers.addRow();
    transitions.push_back(std::move(label)); liveT.push_back(true);
    replacePre(t, std::move(input)); replacePost(t, std::move(output));
    return t;
  }

  /** A transition leaves without a survivor standing for its arcs. */
  void retireTransition(size_t t) {
    if (!liveT[t]) return;
    if (counting) counting->dropArcs("transitions removed whose arcs no survivor stands for");
    dropTransition(t);
  }
  /** `t` duplicates `survivor`, which takes its multiplicity. */
  void fuseTransition(size_t t, size_t survivor) {
    if (!liveT[t]) return;
    if (counting) counting->fused(t, survivor);
    dropTransition(t);
  }
  /** A transition proven never enabled: no arc of the graph was its. */
  void retireDeadTransition(size_t t) {
    if (liveT[t]) dropTransition(t);
  }
  /** The same for a set, one pass per place column they touch. */
  void retireDeadTransitions(std::vector<size_t> ts) {
    std::sort(ts.begin(), ts.end());
    ts.erase(std::unique(ts.begin(), ts.end()), ts.end());
    std::erase_if(ts, [&](size_t t) { return !liveT[t]; });
    if (ts.empty()) return;
    std::vector<size_t> touched;
    for (size_t t : ts)
      for (const auto* mat : {&pre, &post}) {
        const auto& col = mat->getColumn(t);
        for (size_t i = 0; i < col.size(); ++i) touched.push_back(col.keyAt(i));
      }
    std::sort(touched.begin(), touched.end());
    touched.erase(std::unique(touched.begin(), touched.end()), touched.end());
    for (size_t p : touched) { consumers.getColumn(p).removeKeys(ts); producers.getColumn(p).removeKeys(ts); }
    for (size_t t : ts) { pre.getColumn(t).clear(); post.getColumn(t).clear(); liveT[t] = false; ++changes; }
  }
private:
  void dropTransition(size_t t) {
    for (auto pair : {std::pair{&pre, &consumers}, std::pair{&post, &producers}}) {
      auto& col = pair.first->getColumn(t);
      for (size_t i = 0; i < col.size(); ++i) pair.second->getColumn(col.keyAt(i)).put(t, 0);
      col.clear();
    }
    liveT[t] = false; ++changes;
  }
public:

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
    if (counting) throw std::logic_error("Retiring a place under a counting record without saying what it held");
    dropPlace(p);
  }
  /** A constant place leaves, its token total and any free-component
   * coefficient recorded outside the surviving net. */
  void retireConstantPlace(size_t p) {
    if (!liveP[p] || observed[p]) throw std::logic_error("Retiring a protected/inactive place");
    if (counting) counting->constantDropped(p, marks[p]);
    dropPlace(p);
  }
  /** The same for a set, one pass per transition column they touch. */
  void retireConstantPlaces(std::vector<size_t> ps) {
    std::sort(ps.begin(), ps.end());
    ps.erase(std::unique(ps.begin(), ps.end()), ps.end());
    std::erase_if(ps, [&](size_t p) {
      if (!liveP[p]) return true;
      if (observed[p]) throw std::logic_error("Retiring a protected place");
      return false;
    });
    if (ps.empty()) return;
    std::vector<size_t> touched;
    for (size_t p : ps) {
      if (observed[p]) throw std::logic_error("Retiring a protected place");
      for (const auto* mat : {&consumers, &producers}) {
        const auto& col = mat->getColumn(p);
        for (size_t i = 0; i < col.size(); ++i) touched.push_back(col.keyAt(i));
      }
    }
    std::sort(touched.begin(), touched.end());
    touched.erase(std::unique(touched.begin(), touched.end()), touched.end());
    for (size_t t : touched) { pre.getColumn(t).removeKeys(ps); post.getColumn(t).removeKeys(ps); }
    for (size_t p : ps) {
      if (counting) counting->constantDropped(p, marks[p]);
      consumers.getColumn(p).clear(); producers.getColumn(p).clear(); liveP[p] = false; ++changes;
    }
  }
  /** `kept` absorbs `p` (a fused free component); the caller has already
   * summed the markings and moved the arcs. */
  void fusePlace(size_t p, size_t kept) {
    if (counting) counting->placesFused(p, kept);
    dropPlace(p);
  }
private:
  void dropPlace(size_t p) {
    if (!liveP[p] || observed[p]) throw std::logic_error("Retiring a protected/inactive place");
    erasePlaceArcs(p); liveP[p] = false; ++changes;
  }
public:

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
