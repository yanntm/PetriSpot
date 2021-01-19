#pragma once

#include "expr/Expression.h"
#include "expr/Property.h"

#include "ddd/util/ext_hash_map.hh"


namespace d3::util {

template<>
struct hash<petri::expr::Expression *>
{
  size_t
  operator()(const petri::expr::Expression* e) const
  {
    return (size_t)e;
  }
};

template<>
struct equal<petri::expr::Expression *>
{
  bool
  operator()(const petri::expr::Expression* e1,const petri::expr::Expression* e2) const
  {
	  return e1 == e2;
  }
};

} // namespace

namespace petri {
namespace expr {


class AtomicPropManager {
public :
	typedef std::pair<std::string, Expression *> AtomicProp;
	// map place name to index
	typedef ext_hash_map<Expression *, AtomicProp > atommap_t;
private :
	atommap_t atomMap;
	std::vector<AtomicProp> atoms;

	typedef ext_hash_map<std::string, AtomicProp > uniqueMap_t;


	void collectAP(Expression * obj, uniqueMap_t & uniqueMap) {
		if (isPureBool(obj)) {
			// helps to recognize that !AP is the negation of AP
			// Can reduce number of AP as well as help simplifications
			if (obj->getOp() == NOT) {
				obj = obj->childAt(0);
			}
			std::string stringProp = to_string(obj);
			uniqueMap_t::accessor acc;
			bool ins = uniqueMap.insert(acc, stringProp);
			if (ins) {
				acc->second = AtomicProp("p" + std::to_string(atoms.size()), obj);
				atoms.push_back(acc->second);
			}
			atommap_t::accessor acc2;
			atomMap.insert(acc2, obj);
			acc2->second = acc->second;
		} else {
			for (int i=0,ie=obj->nbChildren() ; i < ie ; ++i) {
				collectAP(obj->childAt(i),uniqueMap);
			}
		}
	}


public :

	atommap_t & getAtomMap() {
		return atomMap;
	}

	std::vector<AtomicProp> & getAtoms() {
		return atoms;
	}

	void loadAtomicProps(const std::vector<Property> & props) {
		atoms.clear();
		atomMap.clear();
		uniqueMap_t uniqueMap ;
		// look for atomic propositions
		if (!props.empty()) {
			for (Property prop : props) {
				collectAP(prop.getBody(), uniqueMap);

//				if (prop.getType() == PropertyType.INVARIANT) {
//					Expression be = ((BinOp) prop.getBody()).left;
//					if (prop.getBody().getOp() == Op.EF) {
//						be = Expression.not(be);
//					}
//					atoms.add(new AtomicProp(prop.getName().replaceAll("-", ""), be));
//				} else if (prop.getType() == PropertyType.LTL || prop.getType() == PropertyType.CTL) {
//					collectAP(prop.getBody(), uniqueMap);
//				}
			}

		}

	}

	size_t size() {
		return atoms.size();
	}

	static bool isPureBool(Expression * obj) {
		if (obj == nullptr) {
			return true;
		} else {
			switch (obj->getOp()) {
			case AND:
			case OR:
			case NOT: {
				for (int i = 0, ie = obj->nbChildren(); i < ie; i++) {
					Expression * child = obj->childAt(i);
					if (!isPureBool(child)) {
						return false;
					}
				}
				return true;
			}
			case GT:
			case GEQ:
			case EQ:
			case NEQ:
			case LT:
			case LEQ:
				return true;
			default:
				return false;
			}
		}
	}

	void print(Expression * e, std::ostream & os) {
		atommap_t::const_accessor acc;
		if (atomMap.find(acc, e)) {
			os << acc->second.first ;
		} else {
			switch (e->getOp()) {
			// infix operators
			case AND :
			case OR :
			case U :
			{
				os << "(";
				for (int i=0,ie=e->nbChildren(); i < ie ; i++) {
					if (i > 0) {
						os << " " << to_string(e->getOp()) << " ";
					}
					print(e->childAt(i),os);
				}
				os << ")";
				break;
			}
			// prefix operators
			case G :
			case F :
			case X :
			case NOT :
			{
				os << to_string(e->getOp()) <<"(";
				print(e->childAt(0),os);
				os << ")";
				break;
			}
			default :
				throw "Unexpected operator";
			}
		}
	}

//		ByteArrayOutputStream baos = new ByteArrayOutputStream();
//		{
//			CExpressionPrinter printer = new CExpressionPrinter(new PrintWriter(baos), "src");
//			e.accept(printer);
//			printer.close();
//		}
//		return baos.toString();
//	}
};

}} // namespace
