#pragma once

#include "expr/Expression.h"

namespace petri {
namespace expr {


class VarRef : public Expression {
	int index;

public :
	VarRef(int index) : index(index){
	}

	int eval(const SparseIntArray & state) {
		return state.get(index);
	}

	void print(std::ostream & os)  const  {
		os << "p" << index;
	}
	size_t nbChildren() const {
		return 0;
	}
	Expression * childAt(size_t index) {
		throw "No Children.";
	}

	Op getOp() const {
		return PLACEREF;
	}
	/*
	@Override
	public int getValue() {
		return index;
	}
	
	@Override
	public <T> T accept(ExprVisitor<T> v) {
		return v.visit(this);
	}
	
	@Override
	public String toString() {
		return "s"+index;
	}

	@Override
	public int evalDistance(SparseIntArray state, boolean isNeg) {		
		throw new UnsupportedOperationException();
	}

	@Override
	public Op getOp() {
		return Op.PLACEREF;
	}

	@Override
	public int hashCode() {
		return 2969 * (index +1);
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		VarRef other = (VarRef) obj;
		if (index != other.index)
			return false;
		return true;
	}
	*/
};


}} // namespace
