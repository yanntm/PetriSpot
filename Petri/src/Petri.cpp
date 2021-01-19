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
#include "parse/PTNetLoader.h"
#include "expr/parse/ExprLoader.h"
#include "expr/Expression.h"
#include "expr/AtomicPropManager.h"

#include "Petricube.h"

// Spot's includes, keep it sorted alphabetically
#include <spot/mc/mc_instanciator.hh>
#include <spot/tl/apcollect.hh>
#include <spot/tl/defaultenv.hh>
#include <spot/tl/parse.hh>
#include <spot/twaalgos/dot.hh>
#include <spot/twaalgos/translate.hh>
#include <spot/twacube/twacube.hh>
#include <spot/twacube_algos/convert.hh>

using namespace std;


static bool endsWith(const std::string& str, const std::string& suffix)
{
    return str.size() >= suffix.size() && 0 == str.compare(str.size()-suffix.size(), suffix.size(), suffix);
}


int main(int argc, const char** argv)
{
  if (argc == 1 || argc > 4)
    {
      std::cerr << "usage: petri model.pnml [formula]\n";
      std::cerr << "     model.pnml: the model in the pnml format\n";
      std::cerr << "     [formula]: the LTL formula to verify on model"
		<< "(optional) \n";
      exit(1);
    }


  std::string model(argv[1]);
  std::string formula;
  if (argc >= 3) {
	  formula = argv[2];
  }
  bool ltl = argc == 3;

  bool display = true;

  try {
    SparsePetriNet * pn = loadXML(model);

    if ( endsWith(formula,".xml")) {
    	petri::expr::loadXMLProperties(formula,pn);

    	petri::expr::AtomicPropManager apm;
    	apm.loadAtomicProps(pn->getProperties());

    	std::cout << "Found a total of " << apm.getAtoms().size() << " AP" << std::endl;
    	for (const auto & atom : apm.getAtoms()) {
    		std::cout << atom.first << " : " ;
    		atom.second->print(std::cout);
    		std::cout << std::endl;
    	}

   		// Parse formula...
   		for (auto & property : pn->getProperties()) {
   			// Setup the environment and default dictionnary
   			spot::default_environment& env = spot::default_environment::instance();
   			auto dict = spot::make_bdd_dict();

   			std::ostringstream ostr ;
   			// don't forget to negate formula !
   			ostr << "!(";
   	   		apm.print(property.getBody(), ostr);
   	   		ostr << ")";
   	   		std::cout << "Working on formula " << property.getName() << " :" << ostr.str() << std::endl;
   	   		auto pf = spot::parse_infix_psl(ostr.str(), env, false);
   			auto f = pf.f;
   			// and translate it!
   			spot::translator trans(dict);
   			auto prop = trans.run(&f);

   			// collect atomic propositions
   			spot::atomic_prop_set aps;
   			atomic_prop_collect(f, &aps);


   			std::cout << "\nWorking with the following " <<  aps.size() << " atomic propositions:\n";
   			for (spot::atomic_prop_set::const_iterator ap = aps.begin();
   			     ap != aps.end(); ++ap)
   			  {
   			    std::cout << "   " << ap->ap_name() << '\n';
   			  }

   			// Build the equivalent twacube
   			auto propcube = spot::twa_to_twacube(prop);


   			// FIXME do something with propcube

   			auto* pc = new Petricube(*pn, propcube->ap(),apm);
   			// FIXME do something with petricube
   			std::cout << pc->to_string(pc->initial()) << std::endl;

   			// Instanciate the modelchecking algorithm to use
   		        auto result = spot::ec_instanciator<Petricube_ptr,
   					      SparseIntArray,
   		                              PT_iterator,
   					      SparseIntArray_hash,
   		                              SparseIntArray_equal>
   			  (spot::mc_algorithm::BLOEMEN_EC, pc, propcube, true);


   			std::cout << result.value.at(0) << std::endl;
   			if (result.value.at(0) == spot::mc_rvalue::NOT_EMPTY) {
   				std::cout << "FORMULA " << property.getName() << " FALSE " << " TECHNIQUES TODO" << std::endl ;
   			} else if (result.value.at(0) == spot::mc_rvalue::EMPTY) {
   				std::cout << "FORMULA " << property.getName() << " TRUE " << " TECHNIQUES TODO" << std::endl ;
   			}
   			delete pc;
   		}
   		return 0;
    }

    if (display)
    {
    	std::cout << "PN : " ;
    	pn->getFlowPT().print(std::cout);
    	pn->getFlowTP().print(std::cout);
    	std::cout << std::endl ;
    	for (auto & m : pn->getMarks()) {
    		std::cout << m << ","  ;
    	}
    }
    std::cout << '\n' ;

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


	std::cout << "\nWorking with the following " <<  aps.size() << " atomic propositions:\n";
	for (spot::atomic_prop_set::const_iterator ap = aps.begin();
	     ap != aps.end(); ++ap)
	  {
	    std::cout << "   " << ap->ap_name() << '\n';
	  }

	// Build the equivalent twacube
	auto propcube = spot::twa_to_twacube(prop);


	// FIXME do something with propcube
	petri::expr::AtomicPropManager apm;
	auto* pc = new Petricube(*pn, propcube->ap(),apm);
	// FIXME do something with petricube
	std::cout << pc->to_string(pc->initial()) << std::endl;

	// Instanciate the modelchecking algorithm to use
        auto result = spot::ec_instanciator<Petricube_ptr,
			      SparseIntArray,
                              PT_iterator,
			      SparseIntArray_hash,
                              SparseIntArray_equal>
	  (spot::mc_algorithm::DEADLOCK, pc, propcube, true);

	std::cout << result << std::endl;
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
