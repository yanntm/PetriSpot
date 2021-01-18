
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
	// EF, EG, AF, AG, EX, AX,
	// CTL Binary
	// EU, AU,
	// LTL unary
	F, G, X,
	// LTL binary
	U,  
};

}} // namespace

