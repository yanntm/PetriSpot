/*
 * PTNetLoader.h
 *
 *  Created on: Oct 13, 2020
 *      Author: ythierry
 */

#ifndef EXPRPTNETLOADER_H_
#define EXPRPTNETLOADER_H_

#include <vector>
#include "expr/Expression.h"
#include "SparsePetriNet.h"

namespace petri {
namespace expr {

int loadXMLProperties(std::string filename, SparsePetriNet * spec);


}}



#endif /* PTNETLOADER_H_ */
