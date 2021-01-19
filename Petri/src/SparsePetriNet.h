/*
 * SparsePetriNet.h
 *
 *  Created on: Oct 13, 2020
 *      Author: ythierry
 */

#ifndef SPARSEPETRINET_H_
#define SPARSEPETRINET_H_

#include "MatrixCol.h"
#include "expr/Property.h"


class SparsePetriNet {
	std::string name;
	std::vector<int> marks;
	MatrixCol flowPT;
	MatrixCol flowTP;
	std::vector<std::string> tnames;
	std::vector<std::string> pnames;
	int maxArcValue;
	std::vector<petri::expr::Property *> properties;

	static const int DEBUG = 0;

public :
	SparsePetriNet():name("Petri"),maxArcValue(0) {
	}

	const std::string & getName() const {
		return this->name;
	}

	void setName(const std::string & n) {
		this->name = n;
	}

	int addTransition (const std::string & tname) {
		flowPT.appendColumn(SparseIntArray());
		flowTP.appendColumn(SparseIntArray());
		tnames.emplace_back(tname);
		return tnames.size()-1;
	}

	int addPlace (const std::string & pname, int init) {
		flowPT.addRow();
		flowTP.addRow();
		pnames.emplace_back(pname);
		marks.push_back(init);
		return pnames.size()-1;
	}

	void addPreArc (int p, int t, int val) {
		flowPT.getColumn(t).put(p,val);
		maxArcValue = std::max(maxArcValue, val);
	}

	void addPostArc (int p, int t, int val) {
		flowTP.getColumn(t).put(p,val);
		maxArcValue = std::max(maxArcValue, val);
	}

	int getTransitionCount() const {
		return tnames.size();
	}

	int getPlaceCount() const {
		return pnames.size();
	}

	void setMarking (int pid, int val) {
		marks[pid]=val;
	}

	const std::vector<std::string> & getTnames() const {
		return tnames;
	}

	const std::vector<std::string> & getPnames() const {
		return pnames;
	}

	MatrixCol & getFlowPT() {
		return flowPT;
	}
	MatrixCol & getFlowTP() {
		return flowTP;
	}
	const MatrixCol & getFlowPT() const {
		return flowPT;
	}
	const MatrixCol & getFlowTP() const {
		return flowTP;
	}

	int getMaxArcValue() {
		return maxArcValue;
	}

	const std::vector<int> & getMarks() const {
		return marks;
	}

	std::vector<petri::expr::Property *> & getProperties() {
		return properties;
	}
};


#endif /* SPARSEPETRINET_H_ */
