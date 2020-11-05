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

#include "Petricube.h"

// Spot's includes, keep it sorted alphabetically
#include <spot/tl/apcollect.hh>
#include <spot/tl/defaultenv.hh>
#include <spot/tl/parse.hh>
#include <spot/twaalgos/dot.hh>
#include <spot/twaalgos/translate.hh>
#include <spot/twacube/twacube.hh>
#include <spot/twacube_algos/convert.hh>

using namespace std;

int main(void)
{
  bool display = true;
  bool ltl = true;
  std::string formula = "GF a && GF b"; // FIXME

  try {
    SparsePetriNet * pn = loadXML("model.pnml");

    if (display)
      {
	std::cout << "PN : " ;
	pn->getFlowPT().print(std::cout);
	pn->getFlowTP().print(std::cout);
      }

    if (ltl)
      {
	// Setup the environment and default dictionnary
	spot::default_environment& env = spot::default_environment::instance();
	auto dict = spot::make_bdd_dict();

	// Parse formula...
	auto pf = spot::parse_infix_psl(formula, env, false);
	auto f = pf.f;

	// and translate it!
	spot::translator trans(dict);
	auto prop = trans.run(&f);

	// collect atomic propositions
	spot::atomic_prop_set aps;
	atomic_prop_collect(f, &aps);


	std::cout << "\nWorking with the following atomic propositions:\n";
	for (spot::atomic_prop_set::const_iterator ap = aps.begin();
	     ap != aps.end(); ++ap)
	  {
	    std::cout << "   " << ap->ap_name() << '\n';
	  }

	// Build the equivalent twacube
	auto propcube = spot::twa_to_twacube(prop);

	// FIXME do something with propcube

	auto* pc = new Petricube();
	// FIXME do something with petricube
	delete pc;
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
