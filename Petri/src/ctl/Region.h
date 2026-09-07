/*
 * Region.h
 *
 * The region proof for A[a U b] and A[a W b] from a state: a DFS over the
 * states reachable through states where b fails, each owing a, capped in
 * size; a state failing both, and for U a deadlock or a cycle where b never
 * holds, is a counter-example with the DFS path as evidence. Fires and
 * reverts one cursor in place. See algorithm.md.
 */
#ifndef PETRI_CTL_REGION_H_
#define PETRI_CTL_REGION_H_

#include <cstdint>
#include <unordered_map>
#include <vector>

#include "ctl/Cursor.h"
#include "ctl/Evidence.h"
#include "ctl/Verdict.h"
#include "expr/CtlFormula.h"

namespace petri::ctl
{

template<typename T>
  class Checker;

template<typename T>
  class Region
  {
    using CtlFormula = petri::expr::CtlFormula;

    struct Frame
    {
      std::vector<uint32_t> succ; // the enabled transitions when the state was entered
      size_t next = 0;
      uint32_t via = 0;           // the transition fired from the parent state
    };
    static constexpr uint8_t ON_STACK = 1, DONE = 2;

    /** The state of one run; local to it, since sub-obligations re-enter run. */
    struct Search
    {
      std::vector<Frame> stack;
      std::vector<uint32_t> path;
      std::unordered_map<SparseArray<T>, uint8_t> visited;
    };

    Checker<T> &chk;

    enum class Enter
    {
      Leaf, Expand, Fail, Abort
    };

    /** A state just reached: b closes it, a must hold, the deadlock rule, else it is expanded. */
    Enter enter (Search &sr, Cursor<T> &cur, const CtlFormula &a, const CtlFormula &b, bool until, Evidence &ev,
                 const Budget &budget, uint32_t via)
    {
      Stats &st = chk.stats ();
      auto &visited = sr.visited;
      auto &path = sr.path;
      auto &stack = sr.stack;
      Verdict rb = chk.solve (cur, b, nullptr);
      if (rb == Verdict::Unknown) return Enter::Abort;
      if (rb == Verdict::True) {
        visited[cur.state ()] = DONE;
        ++st.regionStates;
        return Enter::Leaf;
      }
      Evidence kid;
      Verdict ra = chk.solve (cur, a, &kid);
      if (ra == Verdict::Unknown) return Enter::Abort;
      if (ra == Verdict::False) {
        ev.what = "A: path to a state where both sides fail";
        ev.path = path;
        return Enter::Fail;
      }
      if (cur.deadlock ()) {
        if (until) {
          ev.what = "AU: path to a deadlock where the right side never holds";
          ev.path = path;
          return Enter::Fail;
        }
        visited[cur.state ()] = DONE;
        ++st.regionStates;
        return Enter::Leaf;
      }
      if (visited.size () >= budget.regionStates) return Enter::Abort;
      visited[cur.state ()] = ON_STACK;
      ++st.regionStates;
      Frame f;
      f.succ = cur.enabled.transitions ();
      f.via = via;
      stack.push_back (std::move (f));
      return Enter::Expand;
    }

  public:
    explicit Region (Checker<T> &c)
        : chk (c)
    {
    }

    /**
     * Prove A[a U b] (until) or A[a W b] from start. True when the region
     * closes under the budget, False with a counter-example path, Unknown
     * when the budget ends or a sub-obligation stays open. node is the
     * formula the region decides and proved its verdict when the region
     * closes (False when proving the negation of node): every region state
     * is remembered with it, and every state on a counter-example path with
     * the opposite (the suffix is a counter-example there too).
     */
    Verdict run (const Cursor<T> &start, const CtlFormula &a, const CtlFormula &b, bool until, Evidence &ev,
                 const Budget &budget, const CtlFormula &node, Verdict proved)
    {
      ++chk.stats ().regions;
      Cursor<T> cur (start);
      Search sr;
      auto &stack = sr.stack;
      auto &path = sr.path;
      auto &visited = sr.visited;
      auto failed = [&] () {
        Cursor<T> re (start);
        for (uint32_t t : path) {
          chk.remember (re.state (), node, negate (proved));
          re.fire (t);
        }
        return Verdict::False;
      };
      switch (enter (sr, cur, a, b, until, ev, budget, 0)) {
      case Enter::Leaf:
        ev.what = "holds at the state itself";
        return Verdict::True;
      case Enter::Fail: return failed ();
      case Enter::Abort: return Verdict::Unknown;
      case Enter::Expand: break;
      }
      size_t clock = 0;
      while (!stack.empty ()) {
        if ((++clock & 255) == 0 && chk.timedOut ()) return Verdict::Unknown;
        Frame &f = stack.back ();
        if (f.next < f.succ.size ()) {
          uint32_t t = f.succ[f.next++];
          cur.fire (t);
          path.push_back (t);
          auto it = visited.find (cur.state ());
          if (it != visited.end ()) {
            if (it->second == ON_STACK && until) {
              ev.what = "AU: cycle where the right side never holds";
              ev.path = path;
              ev.loopIndex = static_cast<long> (loopStart (cur, path));
              return failed ();
            }
            cur.revert (t);
            path.pop_back ();
            continue;
          }
          switch (enter (sr, cur, a, b, until, ev, budget, t)) {
          case Enter::Fail: return failed ();
          case Enter::Abort: return Verdict::Unknown;
          case Enter::Leaf:
            cur.revert (t);
            path.pop_back ();
            break;
          case Enter::Expand: break;
          }
        } else {
          visited[cur.state ()] = DONE;
          uint32_t via = f.via;
          stack.pop_back ();
          if (!stack.empty ()) {
            cur.revert (via);
            path.pop_back ();
          }
        }
      }
      ev.what = "region of " + std::to_string (visited.size ()) + " states closed";
      for (const auto &kv : visited) chk.remember (kv.first, node, proved);
      return Verdict::True;
    }

  private:
    /** Index on the path of the state the last firing returned to: replay the path from start. */
    size_t loopStart (const Cursor<T> &cur, const std::vector<uint32_t> &path)
    {
      Cursor<T> re (cur);
      for (size_t i = path.size (); i-- > 0;) re.revert (path[i]);
      for (size_t i = 0; i < path.size (); ++i) {
        if (re.state () == cur.state ()) return i;
        re.fire (path[i]);
      }
      return 0;
    }
  };

} // namespace petri::ctl

#endif /* PETRI_CTL_REGION_H_ */
