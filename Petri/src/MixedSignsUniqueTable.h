#ifndef MIXEDSIGNSUNIQUETABLE_H_
#define MIXEDSIGNSUNIQUETABLE_H_

#include <unordered_set>
#include <cstddef>
#include "MatrixCol.h"  // Assumed to define MatrixCol<T> and SparseArray<T>
#include "SparseArray.h" // For SparseArray<T> access

template<typename T>
  class MixedSignsUniqueTable
  {
  public:
    MixedSignsUniqueTable (MatrixCol<T> &matrix)
        : colsB (matrix), table (0, Hash { this }, Equal{ this }),
          attemptedInsertions (0), successfulInsertions (0)
    {
    }

    void clear() {
      table.clear ();
    }

    void erase (size_t index)
    {
      table.erase (index);
    }

    bool insert (size_t index)
    {
      attemptedInsertions++;  // Count every attempt
      auto [it, inserted] = table.insert (index);
      if (inserted) {
        successfulInsertions++;  // Count only successful insertions
      }
      return inserted;
    }

    size_t getAttemptedInsertions () const
    {
      return attemptedInsertions;
    }

    size_t getSuccessfulInsertions () const
    {
      return successfulInsertions;
    }

    size_t size () const
    {
      return table.size ();
    }

  private:
    MatrixCol<T> &colsB;
    size_t hash (size_t index) const
    {
      const SparseArray<T> &col = colsB.getColumn (index);
      size_t h = 17;
      for (size_t i = 0, ie = col.size (); i < ie; ++i) {
        size_t key = col.keyAt (i);
        T val = col.valueAt (i);
        h = h * (val > 0 ? 13 : 31) + key;
      }
      return h;
    }

    bool equal (size_t index1, size_t index2) const
    {
      const SparseArray<T> &a = colsB.getColumn (index1);
      const SparseArray<T> &b = colsB.getColumn (index2);
      if (a.size () != b.size ()) return false;
      for (size_t i = 0, ie = a.size (); i < ie; ++i) {
        if (a.keyAt (i) != b.keyAt (i)
            || ((a.valueAt (i) > 0 ? 1 : -1) * (b.valueAt (i) > 0 ? 1 : -1) < 0)) {
          return false;
        }
      }
      return true;
    }

    struct Hash
    {
      size_t operator() (size_t index) const
      {
        return owner->hash (index);
      }
      const MixedSignsUniqueTable *owner;
    };

    struct Equal
    {
      bool operator() (size_t index1, size_t index2) const
      {
        return owner->equal (index1, index2);
      }
      const MixedSignsUniqueTable *owner;
    };

    std::unordered_set<size_t, Hash, Equal> table;
    size_t attemptedInsertions;    // Total insertion attempts
    size_t successfulInsertions;   // Successful (unique) insertions

  };

#endif // MIXEDSIGNSUNIQUETABLE_H_
