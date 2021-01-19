#pragma once

#include "expr/Expression.h"

namespace petri {
namespace expr {


class BoolConstant : public Expression {
	bool value;
public :
	BoolConstant(bool value) :value(value) {
	}

	int eval(const SparseIntArray & state) {
		return value ? 1 : 0;
	}

	void print(std::ostream & os)  const  {
		if (value)
			os << "true";
		else
			os << "false";
	}

	size_t nbChildren() const {
		return 0;
	}
	Expression * childAt(size_t index) {
		throw "No Children.";
	}

	/*
	@Override
	public <T> T accept(ExprVisitor<T> v) {
		return v.visitBool(this);
	}
	
	@Override
	public String toString() {
		return Boolean.toString(value);
	}

	@Override
	public int evalDistance(SparseIntArray state, boolean isNegated) {		
		if (value ^ isNegated) {
			return 0;
		} else {
			return 1000;
		}
	}
	
	@Override
	public Op getOp() {
		return Op.BOOLCONST;
	}

	@Override
	public int getValue() {
		return eval(null);
	}

	@Override
	public int hashCode() {
		if (value) {
			return 12689;
		} else {
			return 20411;
		}
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		BoolConstant other = (BoolConstant) obj;
		if (value != other.value)
			return false;
		return true;
	}
	*/
	
};

}} // namespace

