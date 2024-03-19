#pragma once


#include "SparsePetriNet.h"
#include <vector>


/**
 * Computes a combined flow matrix, stored with column = transition, while removing any duplicates (e.g. due to test arcs or plain redundancy).
 * Updates tnames that is supposed to initially be empty to set the names of the transitions that were kept.
 * This is so we can reinterpret appropriately the Parikh vectors f so desired.
 * @param sr our Petri net
 * @param tnames empty list that will contain the transition names after call.
 * @param representative the mapping from original transition index to their new representative (many to one/surjection)
 * @return a (reduced, less columns than usual) flow matrix
 */
static MatrixCol computeReducedFlow(const SparsePetriNet & sr, std::vector<int> & tnames, std::vector<int> & representative) {
	MatrixCol sumMatrix (sr.getPlaceCount(), 0);
	{
        
		typedef ext_hash_map<const SparseIntArray *, int> map_t;
		map_t::accessor acc;
        map_t seen;

		for (int i=0 ; i < sr.getFlowPT().getColumnCount() ; i++) {
			SparseIntArray combined = SparseIntArray::sumProd(-1, sr.getFlowPT().getColumn(i), 1, sr.getFlowTP().getColumn(i));

			
			bool found = seen.find(acc,&combined);
			if (!found) {
				seen.insert(acc,&combined);
                sumMatrix.appendColumn(combined);
				tnames.push_back(i);
			}
            representative.push_back((*acc).second);
		}
        // add logger ?
	}
	return sumMatrix;
}
/*
static void printInvariant(Collection<SparseIntArray> invariants, List<String> pnames, List<Integer> initial,
		PrintStream out) {
	for (SparseIntArray rv : invariants) {
		StringBuilder sb = new StringBuilder();
		int sum = printEquation(rv, initial, pnames, sb);
		out.println("inv : " + sb.toString() + " = " + sum);
	}
	out.println("Total of " + invariants.size() + " invariants.");
}

static void printInvariant(Collection<SparseIntArray> invariants, List<String> pnames,
		List<Integer> initial) {
	printInvariant(invariants, pnames, initial, System.out);
}

static int printEquation(SparseIntArray inv, List<Integer> initial, List<String> pnames, StringBuilder sb) {
	boolean first = true;
	int sum = 0;
	for (int i = 0; i < inv.size(); i++) {
		int k = inv.keyAt(i);
		int v = inv.valueAt(i);
		if (v != 0) {
			if (first) {
				if (v < 0) {
					sb.append("-");
					v = -v;
				}
				first = false;
			} else {
				if (v < 0) {
					sb.append(" - ");
					v = -v;
				} else {
					sb.append(" + ");
				}
			}
			if (v != 1) {
				sb.append(v + "*" + pnames.get(k));
			} else {
				sb.append(pnames.get(k));
			}
			sum += initial != null ? v * initial.get(k) : 0;
		}
	}
	return sum;
}
*/
