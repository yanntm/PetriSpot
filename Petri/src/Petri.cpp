/*
  ============================================================================
  Name        : Petri.cpp
  Author      :
  Version     :
  Copyright   : Your copyright notice
  ============================================================================
*/

#include <iostream>
#include "SparseIntArray.h"
#include "MatrixCol.h"
#include "SparsePetriNet.h"
#include "Walker.h"
#include "PTNetLoader.h"

using namespace std;

int main(void)
{
  bool display = true;
  try {
    SparsePetriNet * pn = loadXML("model.pnml");

    if (display)
      {
	std::cout << "PN : " ;
	pn->getFlowPT().print(std::cout);
	pn->getFlowTP().print(std::cout);
      }

    Walker walk (*pn);

    // Use fullrand
    if (walk.runDeadlockDetection(1000000, true, 30))
      {
	std::cout << "Deadlock found !" << std::endl;
	delete pn;
	return 0;
      }
    else
      {
	std::cout << "No deadlock found !" << std::endl;
      }

    // Do not use fullrand
    if (walk.runDeadlockDetection(1000000, false, 30))
      {
	std::cout << "Deadlock found !" << std::endl;
      }
    else
      {
	std::cout << "No deadlock found !" << std::endl;
      }
    delete pn;

  }
  catch (const char* e)
    {
      std::cerr << e << std::endl;
      return 1;
    }

  return 0;
}
