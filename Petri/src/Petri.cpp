
#include <iostream>
#include <string>
#include "Walker.h"
#include "PTNetLoader.h"
#include "InvariantMiddle.h"
#include <vector>
#include <unordered_set>
#include <chrono>

using namespace std;

const string FINDDEADLOCK="--findDeadlock";
const string PFLOW="--Pflows";
const string PSEMIFLOW="--Psemiflows";
const string TFLOW="--Tflows";
const string TSEMIFLOW="--Tsemiflows";
const string PATH="-i";
const string QUIET="-q";


int main(int argc, char * argv[]) {
	std::string logMessage = "Running PetriSpot with arguments : [";
	int i=0;
	for (i = 1; i < argc-1; ++i) {
        	logMessage += std::string(argv[i]) + ", ";
    	}
	logMessage += (argc == 1 ? "Empty" : std::string(argv[i])) + "]";
	InvariantMiddle::writeToLog(logMessage);
	auto runtime = std::chrono::steady_clock::now();
	std::string modelPath;
	bool findDeadlock=false;
	bool pflows=false;
	bool tflows=false;
	bool psemiflows=false;
	bool tsemiflows=false;
	bool invariants=false;
	bool quiet=false;

	if (argc == 1 || argc > 6)
    	{
      		cerr << "usage: petri -i model.pnml [options]\n";
     		exit(1);
    	}
	
	for (int i = 1; i < argc; i++) {
		if (argv[i] == PATH) {
			modelPath = argv[++i];
		} else if (argv[i] == QUIET) {
			quiet = true;
		} else if (argv[i] == FINDDEADLOCK) {
			findDeadlock = true;
		} else if (argv[i] == PFLOW) {
			pflows = true;
			invariants = true;
		} else if (argv[i] == PSEMIFLOW) {
			psemiflows = true;
			invariants = true;
		} else if (argv[i] == TFLOW) {
			tflows = true;
			invariants = true;
		} else if (argv[i] == TSEMIFLOW) {
			tsemiflows = true;
			invariants = true;
		} else {
			std::cout << "[WARNING   ] Option : " << argv[i] << " not recognized" << std::endl;
		}
	}

	try {
		SparsePetriNet * pn = loadXML(modelPath);

// 		std::cout << "PN : " ;
// 		std::cout << "\nPre : " << std::endl;
// 		pn->getFlowPT().print(std::cout);
// 		std::cout << "\nPost : " << std::endl;
// 		pn->getFlowTP().print(std::cout);
// 		std::cout << std::endl;

		if (findDeadlock) 
		{
			Walker walk (*pn);

			if (walk.runDeadlockDetection(1000000, true, 30)) {
				std::cout << "Deadlock found !" << std::endl;
				delete pn;
				return 0;
			} else {
				std::cout << "No deadlock found !" << std::endl;
			}
			if (walk.runDeadlockDetection(1000000, false, 30)) {
				std::cout << "Deadlock found !" << std::endl;
			} else {
				std::cout << "No deadlock found !" << std::endl;
			}	
		}
		
		if (invariants) 
		{
			vector<int> repr;
			MatrixCol sumMatrix = InvariantMiddle::computeReducedFlow(*pn, repr);
			if (pflows || psemiflows) {
				auto time = std::chrono::steady_clock::now();
				unordered_set<SparseIntArray> invar;
				if (pflows) {
					invar = InvariantMiddle::computePInvariants(sumMatrix);
				} else {
					invar = InvariantMiddle::computePInvariants(sumMatrix, true, 120);
				}
				std::cout << "Computed " << invar.size() << " P " << (psemiflows?"semi":"") << "flows in " 
					<< std::chrono::duration_cast<std::chrono::milliseconds>
						(std::chrono::steady_clock::now() - time).count() 
							<< " ms." << std::endl;
				if (!quiet) {
    					InvariantMiddle::printInvariant(invar, pn->getPnames(), (*pn).getMarks());
				}
			}
			if (tflows || tsemiflows) {
				auto time = std::chrono::steady_clock::now();
				unordered_set<SparseIntArray> invarT;
				if (tflows) {
					invarT= InvariantMiddle::computeTinvariants(*pn, sumMatrix, repr,false);	
				} else {
					invarT= InvariantMiddle::computeTinvariants(*pn, sumMatrix, repr,true);	
				}
				std::cout << "Computed " << invarT.size() << " T " << (tsemiflows?"semi":"") << "flows in " 
					<< std::chrono::duration_cast<std::chrono::milliseconds>
						(std::chrono::steady_clock::now() - time).count() 
							<< " ms." << std::endl;
				if (!quiet) {
					std::vector<int> emptyVector;
					InvariantMiddle::printInvariant(invarT, pn->getTnames(), emptyVector);
				}
			}
		}

		delete pn;

	} catch (const char* e) {
		std::cout << e << std::endl;
		return 1;
	}
	
	std::cout << "Total runtime " << std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::steady_clock::now() - runtime).count() << " ms." << std::endl;

	return 0;
}
