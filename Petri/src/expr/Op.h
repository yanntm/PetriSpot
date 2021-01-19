#pragma once

#include <string>

namespace petri {
namespace expr {

enum Op {
	// unary
	NOT,
	// binary logic
	AND, OR,
	// arithmetic
	ADD,
	// these additional operators are unused in MCC context
	DIV, MULT, MOD, MINUS,
	// binary comparisons
	EQ, NEQ, GEQ, GT, LEQ, LT,
	// n-ary property atoms
	ENABLED, CARD, BOUND,
	// 0 arity constants
	CONST, // an integer
	BOOLCONST, // true or false
	DEAD, // a deadlock proposition
	// references to objects of the net
	PLACEREF, // a reference to a place
	TRANSREF, // a reference to a transition of the net 
	// CTL unary
	EF, EG, AF, AG, EX, AX,
	// CTL Binary
	EU, AU,
	// LTL unary
	F, G, X,
	// LTL binary
	U,  
};

inline std::string to_string (const Op & op) {
	switch (op) {
	case NOT :
		return "!";
	case AND :
		return "&&";
	case OR :
		return "||";
	case ADD :
		return "+";
	case DIV :
		return "/";
	case MULT :
		return "*";
	case MOD :
		return "%";
	case MINUS :
		return "-";
	case EQ :
		return "=";
	case NEQ :
		return "!=";
	case GEQ :
		return ">=";
	case GT :
		return ">";
	case LEQ :
		return "<=";
	case LT :
		return "<";
	case ENABLED :
		return "ENABLED";
	case CARD :
		return "CARD";
	case BOUND :
		return "BOUND";
	case CONST :
		return "CONST";
	case BOOLCONST :
		return "BCONST";
	case DEAD :
		return "DEAD";
	case PLACEREF :
		return "PLACE";
	case TRANSREF :
		return "TRANS";
	case EF :
		return "EF";
	case EG :
		return "EG";
	case AF :
		return "AF";
	case AG :
		return "AG";
	case EX :
		return "EX";
	case AX :
		return "AX";
	case EU :
		return "EU";
	case AU :
		return "AU";
	case F :
		return "F";
	case G :
		return "G";
	case X :
		return "X";
	case U :
		return "U";
	default :
		throw "Unknown operator.";
	}
}

}} // namespace

