/*
 * Hint.h
 *
 * Side information attached to a property by its name: a Parikh vector from
 * the state equation (firing counts per transition). Plain data; how a walk
 * uses it is the strategy's business (walk/ParikhStrategy.h).
 */
#ifndef PETRI_EXPR_HINT_H_
#define PETRI_EXPR_HINT_H_

#include <cstddef>
#include <utility>
#include <vector>

namespace petri::expr
{

struct ParikhHint
{
  std::vector<std::pair<size_t, long long>> counts; // transition, firings

  bool empty () const
  {
    return counts.empty ();
  }
  /** Sum of the counts: the length of the suggested firing sequence. */
  long long total () const
  {
    long long s = 0;
    for (const auto &c : counts) s += c.second;
    return s;
  }
};

} // namespace petri::expr

#endif /* PETRI_EXPR_HINT_H_ */
