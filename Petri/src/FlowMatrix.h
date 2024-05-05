#ifndef FLOWMATRIX_H_
#define FLOWMATRIX_H_

#include "MatrixCol.h"

class FlowMatrix {
	// To represent the flow matrix, if we can build it. We use a sparse representation.
	// Map variable index -> Transition index -> update to variable (a relative integer)
private:
	MatrixCol flow ;
	MatrixCol read ;
	MatrixCol flowPT ;
	MatrixCol flowTP ;

public:
	FlowMatrix (int nbvar, int nbtrans) : flow(nbvar, nbtrans), read(nbvar, nbtrans), flowPT(nbvar, nbtrans), flowTP(nbvar, nbtrans) {
	}
	
	void addWriteEffect(int tindex, int vindex, int val) {
		if (val == 0)
			return;
		addToColumn(flow.getColumn(tindex), vindex, val);
		if (val < 0) {
			addToColumn(flowPT.getColumn(tindex), vindex, -val);
		} else {
			addToColumn(flowTP.getColumn(tindex), vindex, val);
		}
	}

private:
	void addToColumn(SparseIntArray column, int vindex, int val) {
		int cur = column.get(vindex);				
		cur+=val;
		column.put(vindex, cur);
	}

public:
	void addReadEffect(int tindex, int vindex, int val) {
		SparseIntArray line = flowPT.getColumn(tindex);
		int cur = line.get(vindex);				
		int max= std::max(cur,val);
		if (max != cur) {
			line.put(vindex, max);
			addToColumn(flowTP.getColumn(tindex), vindex, max - cur);
		}
		read.getColumn(tindex).put(vindex, max);
	}
	
	MatrixCol getIncidenceMatrix() {
		return flow;
	}
	
	MatrixCol getRead() {
		return read;
	}
	
	MatrixCol getFlowPT() {
		return flowPT;
	}
	
	MatrixCol getFlowTP() {
		return flowTP;
	}
	
};

#endif /* FLOWMATRIX_H_ */