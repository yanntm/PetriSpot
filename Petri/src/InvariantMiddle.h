#ifndef INVARIANTMIDDLE_H_
#define INVARIANTMIDDLE_H_


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
            		representative.push_back(acc->second);
		}
        // add logger ?
	}
	return sumMatrix;
}


#endif /* INVARIANTMIDDLE_H_ */