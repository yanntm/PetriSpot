/*
 * Walker.h
 *
 *  Created on: Oct 13, 2020
 *      Author: ythierry
 */

#ifndef WALKER_H_
#define WALKER_H_

#include <random>
#include <thread>

#include "MatrixCol.h"
#include "SparseArray.h"
#include "SparsePetriNet.h"
#include <unordered_map>

template<typename T>
class Walker {
	const SparsePetriNet<T> * sr;
	MatrixCol<T> combFlow;
	MatrixCol<T> tFlowPT;
	int * behaviorMap;
	size_t behaviorCount;
	static const int DEBUG = 1;
public :
	Walker(const SparsePetriNet<T> & ssr) : combFlow(ssr.getPnames().size(),0) {
		this->sr = & ssr;
		typedef std::unordered_map<const SparseArray<T> *, std::vector<size_t>, std::hash<SparseArray<T>*>, std::equal_to<SparseArray<T>*>> map_t;
		map_t effects;

		for (size_t i = 0 ;  i < sr->getFlowPT().getColumnCount() ; i ++) {
			combFlow.appendColumn(SparseArray<T>::sumProd(-1, sr->getFlowPT().getColumn(i), 1, sr->getFlowTP().getColumn(i)));
		}

		for (size_t i = 0 ;  i < sr->getFlowPT().getColumnCount() ; i ++) {
			SparseArray<T> & col = combFlow.getColumn(i);
			typename map_t::iterator it = effects.find(&col);
			if (it != effects.end()) {
				it->second.push_back(i);
			} else {
				std::vector<size_t> r;
				r.push_back(i);
				effects.insert(std::make_pair(&col, r));
			}
		}
		behaviorMap = new int [sr->getTnames().size()];
		int i=0;
		for (const auto & ent : effects) {
			for (const auto & t : ent.second) {
				behaviorMap[t]=i;
			}
			i++;
		}
		behaviorCount = effects.size();
		tFlowPT = sr->getFlowPT().transpose();
	}
	~Walker() {
		delete [] behaviorMap;
	}

private :
	int * computeEnabled(const SparseArray<T> & state) {
		int * list  = new int [sr->getTnames().size()+1];
		memset(list,0, (sr->getTnames().size()+1)* sizeof(int));
		int li = 1;
		for (int t = 0, e =  sr->getTnames().size(); t < e; t++) {
			if (SparseArray<T>::greaterOrEqual(state, sr->getFlowPT().getColumn(t))) {
				list[li++] = t;
			}
		}
		list[0] = li -1 ;
		return list;
	}

	void dropAt (int * enabled, int index) {
		if (index < enabled[0]) {
			enabled[index] = enabled[enabled[0]];
		}
		enabled [0] --;
	}
	void add (int * enabled, int value) {
		enabled[enabled[0]+1] = value;
		enabled [0] ++;
	}

public :
	// we just reached "state" by firing tfired
	void updateEnabled (const SparseArray<T> & state, int * enabled, int tfired) {
		if (combFlow.getColumn(tfired).size() == 0) {
			return ;
		}

		bool * seen = new bool [sr->getTnames().size()];
		memset(seen,0,sr->getTnames().size() * sizeof(bool));
		bool * seenEffects = new bool [behaviorCount];
		memset(seenEffects,0,behaviorCount * sizeof(bool));
		for (int i = enabled[0] ; i  >= 1  ; i--) {
			int t = enabled [i];
			if (seen[t] || seenEffects[behaviorMap[t]]) {
				dropAt(enabled,i);
				continue;
			}

			if (SparseArray<T>::greaterOrEqual(state, sr->getFlowPT().getColumn(t))) {
				seen[t] = true;
				seenEffects[behaviorMap[t]] = true;
				continue;
			} else {
				dropAt(enabled,i);
			}

		}

		// the places fed by this transition
		SparseArray<T> tp = combFlow.getColumn(tfired);
		for (int  pi = 0, pie=tp.size() ; pi < pie ; pi++) {
			size_t p = tp.keyAt(pi);
			if (tp.valueAt(pi) > 0) {
				// the set of transitions taking from this place
				SparseArray<T> col = tFlowPT.getColumn(p);
				for (size_t i = 0 ; i < col.size() ; i++) {
					size_t t = col.keyAt(i);
					if (seen[t] || seenEffects[behaviorMap[t]])
						continue;

					if (combFlow.getColumn(t).size()==0) {
						continue;
					}

					if (SparseArray<T>::greaterOrEqual(state, sr->getFlowPT().getColumn(t))) {
						add(enabled, t);
						seen[t] = true;
						seenEffects[behaviorMap[t]] = true;
					}
				}
			}
		}
		delete [] seen;
		delete [] seenEffects;
	}

	SparseArray<T> fire (int t, const SparseArray<T> & state) {
		// NB no enabling check
		return SparseArray<T>::sumProd(1, state, 1, combFlow.getColumn(t));
	}


	/* Thread-safe function that returns a random number between min and max (exclusive).
	This function takes ~142% the time that calling rand() would take. For this extra
	cost you get a better uniform distribution and thread-safety. */
	static int intRand(const int & min, const int & max) {
	    	static thread_local std::mt19937* generator = nullptr;
	    	if (!generator) generator = new std::mt19937(std::chrono::steady_clock::now().time_since_epoch().count());
	    		std::uniform_int_distribution<int> distribution(min, max-1);
	    	return distribution(*generator);
	}


	bool runDeadlockDetection (long nbSteps, bool fullRand, size_t timeout) {
		using std::chrono::steady_clock;
		auto time = steady_clock::now();
		SparseArray<T> state(sr->getMarks());
		int * list = computeEnabled(state);
		dropEmpty(list);

		int last = -1;

		long nbresets = 0;
		int i=0;
		for ( ; i < nbSteps ; i++) {
			size_t dur = std::chrono::duration_cast<std::chrono::milliseconds>(steady_clock::now() - time).count() + 1;
			if (dur > 1000 * timeout) {
				std::cout << "Interrupted Parikh directed walk after " << i  << "  steps, including " << nbresets<<" resets, run timeout after "<< dur <<" ms. (steps per millisecond=" << (i/dur) << " )" ;
				if (DEBUG >=1) {
					std::cout <<  " reached state ";
					state.print(std::cout);
				}
				delete [] list;
				return false;
			}
			if (list[0] == 0) {
				// includes empty effects
				delete [] list;
				list = computeEnabled(state);
				if (list[0] == 0) {
					std::cout << "Finished Parikh directed walk after " << i << "  steps, including " << nbresets << " resets, run found a deadlock after "<< dur <<" ms. (steps per millisecond=" << (i/dur) <<" )";
					if (DEBUG >=1) {
						std::cout <<  " reached state ";
						state.print(std::cout);
					}
					delete [] list;
					return true;
				} else {
					//System.out.println("Dead end with self loop(s) found at step " + i);
					nbresets ++;
					last = -1;
					state = SparseArray<T>(sr->getMarks());
					delete [] list;
					list = computeEnabled(state);
					dropEmpty(list);
					continue;
				}
			}
			int r = intRand(0, list[0])+1;
			int tfired = list[r];

			bool repeat = shouldRepeatLast(last, state);
			if (repeat) {
				tfired = last;
				// iterate firing
				do {
					state = fire ( tfired, state);
					i++;
				} while (SparseArray<T>::greaterOrEqual(state, sr->getFlowPT().getColumn(tfired)));
				updateEnabled(state, list, tfired);
				last = -1;
				continue;
			}

			if (list[0]==1 || fullRand || intRand(0, 100) >= 60) {
				SparseArray<T> newstate = fire ( tfired, state);
				// NB : discards empty events
				updateEnabled(newstate, list, tfired);
				last = tfired;
				state = newstate;
			} else {
				// heuristically follow a successor with less outgoing edges

				SparseArray<T> * succ = new SparseArray<T>[list[0]];
				for (int ti = 1 ; ti-1 < list[0] ; ti++) {
					succ[ti-1] = fire(list[ti],state);
					i++;
				}
				int minSucc = sr->getTnames().size()+1;
				int mini = -1;
				int * minList = nullptr;
				for (int ti = 0 ; ti < list[0] ; ti++) {
					int* listC = new int [sr->getTnames().size()+1];
					memcpy(listC,list,(list[0]+1)*sizeof(int));
					updateEnabled(succ[ti], listC, list[ti+1]);
					if (listC[0] < minSucc   || (listC[0] == minSucc && intRand(0, 100) >= 50)) {
						minSucc = listC[0];
						mini = ti;
						delete [] minList;
						minList = listC;
					} else {
						delete [] listC;
					}
				}
				state = succ[mini];
				last = list[mini+1];
				delete[] list;
				list = minList;
				delete [] succ;
			}
		}
		size_t dur = std::chrono::duration_cast<std::chrono::milliseconds>(steady_clock::now() - time).count() + 1;// avoid zero divide
		std::cout << "Random " << (fullRand?"":"directed ") << "walk for " << i << "  steps, including " <<nbresets<< " resets, run took "<< dur <<" ms (no deadlock found). (steps per millisecond=" << (i/dur) << " )";
		if (DEBUG >=1) {
			std::cout <<  " reached state ";
			state.print(std::cout);
		}
		delete [] list;
		return false;
	}

	bool shouldRepeatLast(int last, const SparseArray<T> & state) const {
		bool repeat = false;
		if (last != -1 && sr->getFlowPT().getColumn(last).size() > 0 && intRand(0, 100) < 98 && SparseArray<T>::greaterOrEqual(state, sr->getFlowPT().getColumn(last))) {
			// make sure there is no divergent behavior here
			const SparseArray<T> & combb = combFlow.getColumn(last);
			for (int j=0, je = combb.size() ; j < je ; j++) {
				if (combb.valueAt(j) < 0) {
					repeat = true;
					break;
				}
			}
		}
		return repeat;
	}


private:
	/** update a list of enabling to remove empty effect transitions*/
	void dropEmpty(int * enabled) {
		for (int i = enabled[0] ; i  >= 1  ; i--) {
			int t = enabled [i];
			if (combFlow.getColumn(t).size() == 0) {
				dropAt(enabled,i);
			}
		}
	}

	void dropUnavailable (int * enabled, const SparseArray<T> & parikh) {
		for (int i = enabled[0] ; i  >= 1  ; i--) {
			int t = enabled [i];
			if (parikh.get(t) <= 0) {
				dropAt(enabled,i);
			}
		}
	}


};



#endif /* WALKER_H_ */
