#ifndef FLOWMATRIX_H_
#define FLOWMATRIX_H_

#include "MatrixCol.h"

template<typename T>
  class FlowMatrix
  {
    // To represent the flow matrix, if we can build it. We use a sparse representation.
    // Map variable index -> Transition index -> update to variable (a relative integer)
  private:
    MatrixCol<T> flow;
    MatrixCol<T> read;
    MatrixCol<T> flowPT;
    MatrixCol<T> flowTP;

  public:
    FlowMatrix (int nbvar, int nbtrans)
        : flow (nbvar, nbtrans), read (nbvar, nbtrans), flowPT (nbvar, nbtrans), flowTP (
            nbvar, nbtrans)
    {
    }

    void addWriteEffect (int tindex, int vindex, T val)
    {
      if (val == 0) return;
      addToColumn (flow.getColumn (tindex), vindex, val);
      if (val < 0) {
        addToColumn (flowPT.getColumn (tindex), vindex, -val);
      } else {
        addToColumn (flowTP.getColumn (tindex), vindex, val);
      }
    }

  private:
    void addToColumn (SparseArray<T> &column, int vindex, T val)
    {
      T cur = column.get (vindex);
      cur += val;
      column.put (vindex, cur);
    }

  public:
    void addReadEffect (int tindex, int vindex, T val)
    {
      SparseArray<T> &line = flowPT.getColumn (tindex);
      T cur = line.get (vindex);
      T max = std::max (cur, val);
      if (max != cur) {
        line.put (vindex, max);
        addToColumn (flowTP.getColumn (tindex), vindex, max - cur);
      }
      read.getColumn (tindex).put (vindex, max);
    }

    MatrixCol<T> getIncidenceMatrix () const
    {
      return flow;
    }

    MatrixCol<T> getRead () const
    {
      return read;
    }

    MatrixCol<T> getFlowPT () const
    {
      return flowPT;
    }

    MatrixCol<T> getFlowTP () const
    {
      return flowTP;
    }
  };

#endif /* FLOWMATRIX_H_ */
