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
#include "SparseIntArray.h"
#include "MatrixCol.h"
#include "SparsePetriNet.h"
#include "Walker.h"
#include "PTNetLoader.h"
#include "ReducedFlow.h"
#include <vector>

using namespace std;


//	private static final String PFLOW="--Pflows";
//	private static final String PSEMIFLOW="--Psemiflows";
//	private static final String TFLOW="--Tflows";
//	private static final String TSEMIFLOW="--Tsemiflows";


int main(void) {
	bool findDeadlock = false;
	bool pflow=false;
	bool tflow=false;
	bool psemiflow = false;
	bool tsemiflow=false;
	bool invariants=false;

	try {

		// parse -i input.pnml
		const char * modelPath = "model.pnml";

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
		} else if (pflow || tflow || psemiflow || tsemiflow) {
			// Go !
//			if (invariants) {
//				vector<int> tnames;
//				vector<int> repr = new ArrayList<>();
//				MatrixCol sumMatrix = computeReducedFlow(reader.getSPN(), tnames, repr);
//				SparsePetriNet spn = reader.getSPN();
//				if (pflows || psemiflows) {
//					long time = System.currentTimeMillis();
//					Set<SparseIntArray> invar;
//					if (pflows) {
//						invar = InvariantCalculator.computePInvariants(sumMatrix);
//					} else {
//						invar = InvariantCalculator.computePInvariants(sumMatrix, true, 120);
//					}
//					System.out.println("Computed "+invar.size()+" P "+(psemiflows?"semi":"")+" flows in "+(System.currentTimeMillis()-time)+" ms.");
//					InvariantSet inv = new InvariantSet(invar, sumMatrix.transpose());
//					inv.print(System.out, spn.getPnames(), spn.getMarks());
//	//				InvariantCalculator.printInvariant(invar, spn.getPnames(), reader.getSPN().getMarks());
//				}
//
//				if (tflows || tsemiflows) {
//					long time = System.currentTimeMillis();
//					Set<SparseIntArray> invarT;
//					if (tflows) {
//						invarT= DeadlockTester.computeTinvariants(reader.getSPN(), sumMatrix, tnames,false);
//					} else {
//						invarT= DeadlockTester.computeTinvariants(reader.getSPN(), sumMatrix, tnames,true);
//					}
//					List<Integer> empty = new ArrayList<>(tnames.size());
//					for (int i=0 ; i < tnames.size(); i++) empty.add(0);
//					List<String> strtnames = tnames.stream().map(id -> spn.getTnames().get(id)).collect(Collectors.toList());
//					System.out.println("Computed "+invarT.size()+" T "+(psemiflows?"semi":"")+" flows in "+(System.currentTimeMillis()-time)+" ms.");
//					InvariantSet inv = new InvariantSet(invarT, sumMatrix);
//					inv.print(System.out, strtnames, empty);
//					//InvariantCalculator.printInvariant(invarT, strtnames, empty );
//				}
//				SparseIntArray inv = DeadlockTester.findPositiveTsemiflow(sumMatrix);
//
//				return null;
//			}

		}


		delete pn;

	} catch (const char* e) {
		std::cout << e << std::endl;
		return 1;
	}
	return 0;
}
