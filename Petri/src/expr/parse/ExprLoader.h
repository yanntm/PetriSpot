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

namespace petri {
namespace expr {

std::vector<Expression *> loadXML(std::string filename);


}}



#endif /* PTNETLOADER_H_ */
