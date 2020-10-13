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

using namespace std;

int main(void) {
	try {
		SparsePetriNet * pn = loadXML("model.pnml");

//		std::cout << "PN : " ;
//		pn->getFlowPT().print(std::cout);
//		pn->getFlowTP().print(std::cout);

		Walker walk (*pn);

		if (walk.runDeadlockDetection(1000000, true, 30)) {
			std::cout << "Deadlock found !" << std::endl;
			return 0;
		} else {
			std::cout << "No deadlock found !" << std::endl;
		}
		if (walk.runDeadlockDetection(1000000, false, 30)) {
			std::cout << "Deadlock found !" << std::endl;
		} else {
			std::cout << "No deadlock found !" << std::endl;
		}
		delete pn;

	} catch (const char* e) {
		std::cout << e << std::endl;
		return 1;
	}
	return 0;
}
