/*
 * Evidence.h
 *
 * The witness tree of a verdict: at each node what was shown, the firing
 * sequence from the node's state, a loop index for a lasso, and the
 * sub-evidence for the children the path relied on.
 */
#ifndef PETRI_CTL_EVIDENCE_H_
#define PETRI_CTL_EVIDENCE_H_

#include <cstdint>
#include <ostream>
#include <string>
#include <vector>

namespace petri::ctl
{

struct Evidence
{
  std::string what;           // "EF: path to b", "AX: successor fails", "region of 12 states"...
  std::vector<uint32_t> path; // transitions fired from the node's state, in order
  long loopIndex = -1;        // lasso: the path's state (0 = start) the last firing returns to
  std::vector<Evidence> kids;

  static Evidence of (std::string what)
  {
    Evidence e;
    e.what = std::move (what);
    return e;
  }

  void print (std::ostream &os, const std::vector<std::string> *tnames, size_t indent = 0) const
  {
    os << std::string (indent, ' ') << what;
    if (!path.empty ()) {
      os << " :";
      for (size_t i = 0; i < path.size (); ++i) {
        if (loopIndex >= 0 && static_cast<size_t> (loopIndex) == i) os << " [";
        os << " " << (tnames ? (*tnames)[path[i]] : "t" + std::to_string (path[i]));
      }
      if (loopIndex >= 0) os << " ]*";
    }
    os << "\n";
    for (const auto &k : kids) k.print (os, tnames, indent + 2);
  }
};

} // namespace petri::ctl

#endif /* PETRI_CTL_EVIDENCE_H_ */
