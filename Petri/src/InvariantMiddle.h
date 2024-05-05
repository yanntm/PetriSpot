#ifndef INVARIANTMIDDLE_H_
#define INVARIANTMIDDLE_H_


#include "SparsePetriNet.h"
#include "InvariantCalculator.h"
#include <chrono>
#include <thread>
#include <future>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>

/**
 * A front-end for functionality computing invariants. Underlying code is
 * adapted from CvO-Theory group's APT : https://github.com/CvO-Theory/apt 
 * See also uniol.apt package and classes.
 * 
 * @author ythierry
 *
 */
class InvariantMiddle {
private:
	static MatrixCol last;
	static std::mutex lock;
	static std::unordered_set<SparseIntArray> lastInv;
	
public:
    	static const int DEBUG = 0;

    	/**
    	 * Guaranteed polynomial runtime, returns flows (with positive AND negative
    	 * coefficients)
    	 * 
    	 * @param pn representing the Petri net approximation
   	  * @return a set of invariants, i.e. coeffs for each variable such that the sum
    	 *         is constant in all markings/states.
    	 */
    	static std::unordered_set<SparseIntArray> computePInvariants(FlowMatrix pn) {
        	std::unordered_set<SparseIntArray> invar;
        	auto time = std::chrono::steady_clock::now();

		std::ofstream logFile("log.txt", std::ios::app);
		if (!logFile.is_open()) {
        		std::cerr << "Erreur : Impossible d'ouvrir le fichier de journalisation." << std::endl;
    		}

        	try {
            		invar = InvariantCalculator::calcSInvariants(pn, InvariantCalculator::InvariantAlgorithm::PIPE, false);
            		logFile << "Computed " << invar.size() << " place invariants in "
                      		<< std::chrono::duration_cast<std::chrono::milliseconds>
					(std::chrono::steady_clock::now() - time).count()
                      				<< " ms" << std::endl;
        	} catch (std::overflow_error& e) {
            		invar = std::unordered_set<SparseIntArray>();
            		logFile << "Invariants computation overflowed in "
                      		<< std::chrono::duration_cast<std::chrono::milliseconds>
					(std::chrono::steady_clock::now() - time).count()
                      				<< " ms" << std::endl;
        	}
		logFile.close();
        	return invar;
    	}

	static void printInvariant(const std::unordered_set<SparseIntArray>& invariants, const std::vector<std::string>& pnames, const std::vector<int>& initial, std::ostream& out) {
        	for (const auto& rv : invariants) {
            		std::stringstream sb;
            		try {
                		long sum = printEquation(rv, initial, pnames, sb);          
                		out << "inv : " << sb.str() << " = " << sum << std::endl;
            		} catch (std::overflow_error& e) {
                		std::cerr << "Overflow of 'long' when computing constant for invariant." << std::endl;
            		}
        	}
        	out << "Total of " << invariants.size() << " invariants." << std::endl;
    	}

	static void printInvariant(const std::unordered_set<SparseIntArray>& invariants, const std::vector<std::string>& pnames, const std::vector<int>& initial) {
        	printInvariant(invariants, pnames, initial, std::cout);
    	}

	static long printEquation(const SparseIntArray& inv, const std::vector<int>& initial, const std::vector<std::string>& pnames, std::stringstream& sb) {
        	bool first = true;
        	long sum = 0;
        	for (size_t i = 0; i < inv.size(); i++) {
            		int k = inv.keyAt(i);
            		int v = inv.valueAt(i); 
            		if (v != 0) {
                		if (!first) {
                    			if (v < 0) {
                        			sb << " - ";
                        			v = -v;
                    			} else {
                        			sb << " + ";
                    			}
                		} else {
                    			if (v < 0) {
                        			sb << "-";
                        			v = -v;
                    			}
                    			first = false;
                		}
                		if (v != 1) {
                    			sb << v << "*" << pnames[k];
                		} else {
                    			sb << pnames[k];
                		}
                		if (!initial.empty()) {
                    			sum = addExact(sum, (multiplyExact((long)v, initial[k])));	
                		}
            		}
        	}
	        return sum;
    	}

private:
	static int32_t addExact(int32_t x, int32_t y) {
		int32_t result;
		if (__builtin_add_overflow(x, y, &result)) {
			throw std::overflow_error("Overflow in addition");
		}
		return result;
	}

	static int32_t multiplyExact(int32_t x, int32_t y) {
		int32_t result;
		if (__builtin_mul_overflow(x, y, &result)) {
			throw std::overflow_error("Overflow in multiplication");
		}
		return result;
	}

public:
	static std::unordered_set<SparseIntArray> computePInvariants(const MatrixCol& pn) {
        	return computePInvariants(pn, false, 120);
    	}

    	static std::unordered_set<SparseIntArray> computePInvariants(const MatrixCol& pn, bool onlyPositive, int timeout) {
        	std::vector<std::thread> threads;
        	std::promise<std::unordered_set<SparseIntArray>> promise;
        	auto future = promise.get_future();

        	threads.emplace_back([&](){
            		std::unordered_set<SparseIntArray> result = computePInvariants(pn, onlyPositive);
            		promise.set_value(result);
        	});

        	auto status = future.wait_for(std::chrono::seconds(timeout));
        	if (status == std::future_status::timeout) {
			std::ofstream logFile("log.txt", std::ios::app);
			if (!logFile.is_open()) {
        			std::cerr << "Erreur : Impossible d'ouvrir le fichier de journalisation." << std::endl;
    			}
            		for (auto& thread : threads) {
                		if (thread.joinable()) {
                    			thread.detach();
                		}
            		}
            		logFile << "Invariant computation timed out after " << timeout << " seconds." << std::endl;
			logFile.close();
            		return std::unordered_set<SparseIntArray>();
        	} else if (status == std::future_status::ready) {
            		return future.get();
        	} else {
			std::ofstream logFile("log.txt", std::ios::app);
			if (!logFile.is_open()) {
        			std::cerr << "Erreur : Impossible d'ouvrir le fichier de journalisation." << std::endl;
    			}
            		for (auto& thread : threads) {
                		if (thread.joinable()) {
                    			thread.detach();
                		}
            		}
            		logFile << "Error: Future is in invalid state." << std::endl;
			logFile.close();
            		return std::unordered_set<SparseIntArray>();
        	}
    	}

private:
	static void cache(const MatrixCol& pn, const std::unordered_set<SparseIntArray>& inv) {
		std::lock_guard<std::mutex> guard(lock);
		last = pn;
		lastInv = inv;
	}

	static std::unordered_set<SparseIntArray> checkCache(const MatrixCol& pn) {
        	std::lock_guard<std::mutex> guard(lock);
		if (pn.equals(last)) {
			std::ofstream logFile("log.txt", std::ios::app);
			if (!logFile.is_open()) {
        			std::cerr << "Erreur : Impossible d'ouvrir le fichier de journalisation." << std::endl;
    			}
			logFile << "Invariant cache hit." << std::endl;
            		logFile.close();
			return lastInv;
		} else {
			return std::unordered_set<SparseIntArray>();
		}
	}

public:
	static std::unordered_set<SparseIntArray> computeTinvariants(const SparsePetriNet & sr, const MatrixCol & sumMatrix, const std::vector<int>& repr,
			bool onlyPositive) {

		std::unordered_map<int, std::vector<int>> repSet = computeMap(repr);
		std::unordered_set<SparseIntArray> invarT = computePInvariants(sumMatrix.transpose(), onlyPositive);

		if (DEBUG >= 1 && !invarT.empty()) {
			std::vector<int> empty(sr.getTransitionCount(), 0);
        		printInvariant(invarT, sr.getTnames(), empty);
    		}
		// so we have T invariants, using the reduced flow matrix
		std::unordered_set<SparseIntArray> reindexT;
		// reinterpret over the original indexes of transitions
		for (const auto& inv : invarT) {
        		std::vector<SparseIntArray> toadd = {SparseIntArray()};
        		for (size_t i = 0; i < inv.size(); ++i) {
            			int t = inv.keyAt(i);
            			int val = inv.valueAt(i);
            			const std::vector<int>& images = repSet[t];
            			if (images.size() > 1) {
                			for (const auto& img : images) {
               	 				std::vector<SparseIntArray> toadd2;
                    				for (const auto& b : toadd) {
                        				SparseIntArray mod = b;
                        				mod.put(img, val);
                        				toadd2.push_back(mod);
                    				}
						toadd = std::move(toadd2);
                			}
            			} else {
                			for (auto& b : toadd) {
                    				b.put(images[0], val);
                			}
            			}
        		}
        		for (const auto& b : toadd) {
            			reindexT.insert(b);
        		}
    		}

		if (DEBUG >= 2 && !invarT.empty()) {
        		printInvariant(reindexT, sr.getTnames(), {});
    		}

		return invarT;
	}

	static std::unordered_set<SparseIntArray> computeTinvariants(const SparsePetriNet & sr, MatrixCol sumMatrix,
			const std::vector<int> & repr) {
		return computeTinvariants(sr, sumMatrix, repr, true);
	}

	static std::unordered_set<SparseIntArray> computePInvariants(const MatrixCol& pn, bool onlyPositive) {
    		std::unordered_set<SparseIntArray> invar = checkCache(pn);
    		if (!invar.empty()) {
        		return invar;
    		}

    		auto startTime = std::chrono::steady_clock::now();
    		try {
        		invar = InvariantCalculator::calcInvariantsPIPE(pn.transpose(), onlyPositive);
        		cache(pn, invar);
			std::ofstream logFile("log.txt", std::ios::app);
			if (!logFile.is_open()) {
        			std::cerr << "Erreur : Impossible d'ouvrir le fichier de journalisation." << std::endl;
    			}
       			logFile << "Computed " << std::to_string(invar.size()) << " invariants in " 
				<< std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
					(std::chrono::steady_clock::now() - startTime).count()) 
						<< " ms" << std::endl;
			logFile.close();
    		} catch (std::exception& e) {
        		invar.clear();
			std::ofstream logFile("log.txt", std::ios::app);
			if (!logFile.is_open()) {
        			std::cerr << "Erreur : Impossible d'ouvrir le fichier de journalisation." << std::endl;
    			}
        		logFile << "Invariants computation overflowed in " 
				<< std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
					(std::chrono::steady_clock::now() - startTime).count()) 
						<< " ms" << std::endl;
			logFile.close();
		}
    		return invar;
	}

	/**
	 * Worst case exponential (time and memory), returns semi-flows (with positive
	 * coefficients only) which are reputed easier to interpret.
	 * 
	 * @param pn representing the Petri net approximation
	 * @return a set of invariants, i.e. coeffs for each variable such that the sum
	 *         is constant in all markings/states.
	 */
	static std::unordered_set<SparseIntArray> computePSemiFlows(FlowMatrix pn) {
		return InvariantCalculator::calcSInvariants(pn, InvariantCalculator::InvariantAlgorithm::PIPE, true);
	}

	/**
	 * Computes a combined flow matrix, stored with column = transition, while
	 * removing any duplicates (e.g. due to test arcs or plain redundancy). Updates
	 * tnames that is supposed to initially be empty to set the names of the
	 * transitions that were kept. This is so we can reinterpret appropriately the
	 * Parikh vectors f so desired.
	 * 
	 * @param sr             our Petri net
	 * @param representative the mapping from original transition index to their new
	 *                       representative (many to one/surjection)
	 * @return a (reduced, less columns than usual) flow matrix
	 */
	static MatrixCol computeReducedFlow(const SparsePetriNet & sr, std::vector<int> & representative) {
		MatrixCol sumMatrix (sr.getPlaceCount(), 0);
		sumMatrix.reserveColumns(sr.getTransitionCount());
		{
			typedef std::unordered_map<const SparseIntArray *, int, std::hash<SparseIntArray*>, std::equal_to<SparseIntArray*>> map_t;
			map_t seen;
			int discarded = 0;
			int curr = 0;
			
			for (int i=0 ; i < sr.getTransitionCount() ; i++) {
				// effects of t
				SparseIntArray combined = SparseIntArray::sumProd(-1, sr.getFlowPT().getColumn(i), 1, sr.getFlowTP().getColumn(i));
				// have we seen this effect ?
				map_t::iterator it = seen.find(&combined);
				if (it == seen.end()) {
					// a new effect
					sumMatrix.appendColumn(combined);
					seen.insert(std::make_pair(& sumMatrix.getColumn(sumMatrix.getColumnCount()-1), curr));
					// this transition is its own representative
					representative.push_back(curr);
					curr++;
				} else {
					// this transition is represented by the element at index :
					representative.push_back(it->second);
					discarded++;
				}
			}
			if (discarded > 0) {
				std::ofstream logFile("log.txt", std::ios::app);
				if (!logFile.is_open()) {
        				std::cerr << "Erreur : Impossible d'ouvrir le fichier de journalisation." << std::endl;
    				}
				logFile << "Flow matrix only has " << sumMatrix.getColumnCount() << " transitions (discarded " << discarded << " similar events)" << std::endl;
				logFile.close();
			}
		}
		return sumMatrix;
	}

	static SparseIntArray transformParikh(const SparseIntArray & parikhori, const std::unordered_map<int, std::vector<int>>& repr) {
		SparseIntArray parikh;
		for (int i = 0, e = parikhori.size(); i < e; i++) {
			int t = parikhori.keyAt(i);
			int k = parikhori.valueAt(i);
			auto it = repr.find(t);
        		if (it != repr.end()) {
            			for (int tr : it->second) {
                			parikh.put(tr, k);
            			}
        		}
		}
		return parikh;
	}

	static std::unordered_map<int, std::vector<int>> computeMap(const std::vector<int>& repr) {
    		std::unordered_map<int, std::vector<int>> repSet;
    		for (size_t i = 0; i < repr.size(); ++i) {
        		int t = i;
        		repSet[repr[t]].push_back(t);
    		}
    		return repSet;
	}

};

#endif /* INVARIANTMIDDLE_H_ */
