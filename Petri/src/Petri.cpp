/*
 ============================================================================
 Name        : Petri.cpp
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C++,
 ============================================================================
 */

#include <iostream>
#include <string.h>
#include "SparseIntArray.h"
#include "MatrixCol.h"
#include "SparsePetriNet.h"
#include "Walker.h"
#include "PTNetLoader.h"
#include "ReducedFlow.h"
#include <vector>
#include <set>

using namespace std;

const string FINDDEADLOCK="--findDeadlock";
const string PFLOWS="--Pflows";
const string PSEMIFLOWS="--Psemiflows";
const string TFLOWS="--Tflows";
const string TSEMIFLOWS="--Tsemiflows";

int main(int argc, char * argv[]) {
	bool findDeadlock=false;
	bool pflows=false;
	bool tflows=false;
	bool psemiflows=false;
	bool tsemiflows=false;
	bool invariants=false;

	if (argc == 1 || argc > 4)
    	{
      	std::cerr << "usage: petri model.pnml [flags]\n";
      	std::cerr << "     model.pnml: the model in the pnml format\n";
      	std::cerr << "     [flags]: findDeadlock or flows to compute"
			<< " (optional) \n";
     	exit(1);
    	}
	
	std::string modelPath(argv[1]);
	if (argc > 2) {
		if (argv[2] == FINDDEADLOCK) {
			findDeadlock = true;
		} else if (argv[2] == PFLOWS || argv[2] == PSEMIFLOWS || argv[2] == TFLOWS || argv[2] == TSEMIFLOWS) {
			invariants = true;
			if (argv[2] == PFLOWS) {
				pflows = true;
			} else if (argv[2] == PSEMIFLOWS) {
				psemiflows = true;
			} else if (argv[2] == TFLOWS) {
				tflows = true;
			} else if (argv[2] == TSEMIFLOWS) {
				tsemiflows = true;
			}
		}
	}
	if (argc > 3) {
		if (argv[3] == TFLOWS) {
			tflows = true;
		} else if (argv[3] == TSEMIFLOWS) {
			tsemiflows = true;
		}
	}
		

	try {
		
		SparsePetriNet * pn = loadXML(modelPath);

		std::cout << "PN : " ;
		std::cout << "\nPre : " << std::endl;
		pn->getFlowPT().print(std::cout);
		std::cout << "\nPost : " << std::endl;
		pn->getFlowTP().print(std::cout);
		std::cout << std::endl;

		if (findDeadlock) {

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
		} else if (pflows || tflows || psemiflows || tsemiflows) {
			// Go !
			if (invariants) {
				vector<int> tnames;
				vector<int> repr;
				MatrixCol sumMatrix = computeReducedFlow(*pn, tnames, repr);
				if (pflows || psemiflows) {
//					long time = System.currentTimeMillis();
					set<SparseIntArray> invar;
					if (pflows) {
//						invar = InvariantCalculator::calcInvariantsPIPE(sumMatrix);
					} else {
//						invar = computePInvariants(sumMatrix, true, 120);
					}
//					std::cout << "Computed " << invar.size() << " P " << (psemiflows?"semi":"") << " flows in " << " ms." << endl;
//					InvariantSet inv = new InvariantSet(invar, sumMatrix.transpose());
//					inv.print(System.out, spn.getPnames(), spn.getMarks());
//    				printInvariant(invar, pn->getPnames(), (*pn).getMarks());
				}
//
				if (tflows || tsemiflows) {
//					long time = System.currentTimeMillis();
					set<SparseIntArray> invarT;
					if (tflows) {
//						invarT= DeadlockTester.computeTinvariants(reader.getSPN(), sumMatrix, tnames,false);
					} else {
//						invarT= DeadlockTester.computeTinvariants(reader.getSPN(), sumMatrix, tnames,true);
					}
					vector<int> empty(tnames.size());
//					for (int i=0 ; i < tnames.size(); i++) empty.add(0);
//					vector<string> strtnames = tnames.stream().map(id -> spn.getTnames().get(id)).collect(Collectors.toList());
//					System.out.println("Computed "+invarT.size()+" T "+(psemiflows?"semi":"")+" flows in "+(System.currentTimeMillis()-time)+" ms.");
//					InvariantSet inv = new InvariantSet(invarT, sumMatrix);
//					inv.print(System.out, strtnames, empty);
//					//InvariantCalculator.printInvariant(invarT, strtnames, empty );
				}
//				SparseIntArray inv = DeadlockTester.findPositiveTsemiflow(sumMatrix);
//				
//				delete pn;
//				return null;
			}

		}


		delete pn;

	} catch (const char* e) {
		std::cout << e << std::endl;
		return 1;
	}
	return 0;
}
