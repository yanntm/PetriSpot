#ifndef INVARIANTMIDDLE_H_
#define INVARIANTMIDDLE_H_


#include "SparsePetriNet.h"
#include <vector>
#include <unordered_map>

/**
 * Computes a combined flow matrix, stored with column = transition, while removing any duplicates (e.g. due to test arcs or plain redundancy).
 * Updates tnames that is supposed to initially be empty to set the names of the transitions that were kept.
 * This is so we can reinterpret appropriately the Parikh vectors f so desired.
 * @param sr our Petri net
 * @param tnames empty list that will contain the transition names after call.
 * @param representative the mapping from original transition index to their new representative (many to one/surjection)
 * @return a (reduced, less columns than usual) flow matrix
 */
static MatrixCol computeReducedFlow(const SparsePetriNet & sr, std::vector<int> & representative) {
	MatrixCol sumMatrix (sr.getPlaceCount(), 0);
	sumMatrix.reserveColumns(sr.getTransitionCount());
	{

		typedef std::unordered_map<const SparseIntArray *, int, std::hash<SparseIntArray*>, std::equal_to<SparseIntArray*>> map_t;

		map_t seen;

		for (int i=0 ; i < sr.getTransitionCount() ; i++) {
			// effects of t
			SparseIntArray combined = SparseIntArray::sumProd(-1, sr.getFlowPT().getColumn(i), 1, sr.getFlowTP().getColumn(i));

			std::cout << "Transition " << i << " effect " << combined << std::endl;
			// have we seen this effect ?
			map_t::iterator it = seen.find(&combined);

			if (it == seen.end()) {
				// a new effect

				sumMatrix.appendColumn(combined);

				std::cout << "Transition " << i << " after append " << sumMatrix.getColumn(sumMatrix.getColumnCount()-1) << std::endl;

				seen.insert(std::make_pair(& sumMatrix.getColumn(sumMatrix.getColumnCount()-1) ,i));

				// this transition is its own representative
				representative.push_back(i);
			} else {
				// this transition is represented by the element at index :
				representative.push_back(it->second);
			}
		}
		auto kv1 = * seen.begin();
		std::cout <<"REFERENCE :" << "key =" << kv1.first << " v=" << kv1.second << " hash=" << std::hash<SparseIntArray*>()(kv1.first) << " effect " << *kv1.first << std::endl;
		for (auto & kv : seen) {
			std::cout << "key =" << kv.first << " v=" << kv.second << " hash=" << std::hash<SparseIntArray*>()(kv.first) << " effect " << *kv.first << std::endl;
			std::cout << "Are equal to reference ?" << (*kv1.first == *kv.first) << " std::equal ?" << std::equal_to<SparseIntArray*>()(kv1.first,kv.first) << std::endl;
		}

        // add logger ?
		size_t b4 = sr.getTransitionCount();
		size_t after = sumMatrix.getColumnCount();
		std::cout << "Initially " << b4 << " transitions. After " << after << " remain. Removed " << (b4-after) << "transitions having same effects." << std::endl;
	}
	return sumMatrix;
}


#endif /* INVARIANTMIDDLE_H_ */
