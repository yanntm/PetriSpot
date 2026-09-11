#pragma once
#include "reduction/Workspace.h"
#include "core/Arithmetic.hpp"
namespace petri::reduction {
template<class T> bool covers(const SparseArray<T>& a, const SparseArray<T>& b) {
  for (size_t i = 0; i < b.size(); ++i) if (a.get(b.keyAt(i)) < b.valueAt(i)) return false;
  return true;
}
template<class T> SparseArray<T> transitionEffect(const Workspace<T>& w, size_t t) {
  SparseArray<T> result;
  const auto& pre = w.pre.getColumn(t); const auto& post = w.post.getColumn(t);
  for (size_t i = 0; i < pre.size(); ++i) result.put(pre.keyAt(i), -pre.valueAt(i));
  for (size_t i = 0; i < post.size(); ++i)
    result.put(post.keyAt(i), petri::addExact(result.get(post.keyAt(i)), post.valueAt(i)));
  return result;
}
template<class T> SparseArray<T> addColumns(SparseArray<T> a, const SparseArray<T>& b) {
  for (size_t i = 0; i < b.size(); ++i) a.put(b.keyAt(i), petri::addExact(a.get(b.keyAt(i)), b.valueAt(i)));
  return a;
}
}
