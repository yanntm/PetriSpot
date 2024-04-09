
#include <iostream>
#include <string>
#include "SparseIntArray.h"
#include "MatrixCol.h"
#include "SparsePetriNet.h"
// #include "Walker.h"
#include "PTNetLoader.h"
#include "InvariantMiddle.h"
#include "InvariantCalculator.h"
#include <vector>
#include <unordered_set>

using namespace std;

const string FINDDEADLOCK="--findDeadlock";
const string PFLOW="--Pflows";
const string PSEMIFLOW="--Psemiflows";
const string TFLOW="--Tflows";
const string TSEMIFLOW="--Tsemiflows";
const string INVARIANT="--invariant";


int main(int argc, char * argv[]) {
	bool findDeadlock=false;
	bool pflows=false;
	bool tflows=false;
	bool psemiflows=false;
	bool tsemiflows=false;
	bool invariants=false;
	std::string modelPath;

	if (argc == 1 || argc > 4)
    	{
      	std::cerr << "usage: petri -i model.pnml [flags]\n";
      	std::cerr << "     model.pnml: the model in the pnml format\n";
      	std::cerr << "     [flags]: findDeadlock or flows to compute"
			<< " (optional) \n";
     	exit(1);
    	}
	
	for (int i = 1; i < argc; i++) {
		cout << argv[i] << endl;
		if (strcmp(argv[i], "-i") == 0) {
			modelPath = argv[++i];
		}
		else if (argv[i] == FINDDEADLOCK) {
			findDeadlock = true;
		}
		else if (argv[i] == INVARIANT) {
			invariants = true;
		}
		else if (argv[i] == PFLOW) {
			pflows = true;
			invariants = true;
		}
		else if (argv[i] == PSEMIFLOW) {
			psemiflows = true;
			invariants = true;
		}
		else if (argv[i] == TFLOW) {
			tflows = true;
			invariants = true;
		}
		else if (argv[i] == TSEMIFLOW) {
			tsemiflows = true;
			invariants = true;
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

//			Walker walk (*pn);
//
//			if (walk.runDeadlockDetection(1000000, true, 30)) {
//				std::cout << "Deadlock found !" << std::endl;
//				delete pn;
//				return 0;
//			} else {
//				std::cout << "No deadlock found !" << std::endl;
//			}
//			if (walk.runDeadlockDetection(1000000, false, 30)) {
//				std::cout << "Deadlock found !" << std::endl;
//			} else {
//				std::cout << "No deadlock found !" << std::endl;
//			}
//
//			delete pn;
			return 0;
			
		}
		
		if (invariants) {
			vector<int> repr;
			MatrixCol sumMatrix = computeReducedFlow(*pn, repr);
			if (pflows || psemiflows) {
//				long time = System.currentTimeMillis();
				unordered_set<SparseIntArray> invar;
				if (pflows) {
					invar = InvariantCalculator::calcInvariantsPIPE(sumMatrix.transpose(), false);
				} else {
					invar = InvariantCalculator::calcInvariantsPIPE(sumMatrix.transpose(), true);
				}
				std::cout << "Computed " << invar.size() << " P " << (psemiflows?"semi":"") << " flows in " << " ms." << endl;
////			InvariantSet inv = new InvariantSet(invar, sumMatrix.transpose());
////			inv.print(System.out, spn.getPnames(), spn.getMarks());
//    			printInvariant(invar, pn->getPnames(), (*pn).getMarks());
			}
			if (tflows || tsemiflows) {
//				long time = System.currentTimeMillis();
				unordered_set<SparseIntArray> invarT;
				if (tflows) {
//					invarT = DeadlockTester.computeTinvariants(reader.getSPN(), sumMatrix, repr,false);
				} else {
//					invarT = DeadlockTester.computeTinvariants(reader.getSPN(), sumMatrix, repr,true);
				}
				cout << "Computed " << invarT.size() << " T " << (psemiflows?"semi":"") << " flows in " << " ms." << endl;
////			InvariantSet inv = new InvariantSet(invarT, sumMatrix);
////			inv.print(System.out, strtnames, empty);
//				InvariantCalculator.printInvariant(invarT, strtnames, empty );
			}
////		SparseIntArray inv = DeadlockTester.findPositiveTsemiflow(sumMatrix);
		}

		delete pn;

	} catch (const char* e) {
		std::cout << e << std::endl;
		return 1;
	}
	return 0;
}
