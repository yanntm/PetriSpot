#ifndef PTNETHANDLER_H_
#define PTNETHANDLER_H_

#include "SparsePetriNet.h"
#include "expr/NaryOp.h"
#include "expr/BinOp.h"
#include "expr/VarRef.h"
#include "expr/Constant.h"
#include "ddd/util/ext_hash_map.hh"

#include <expat.h>
#include <stack>


namespace petri::expr {

	/**
	 * A class to parse a PT model from an PNML file.
	 * @author Yann Thierry-Mieg 2015
	 */
	class ExprHandler {

		// context stack
		std::stack<void *> stack;

		// object constructed
		std::vector<Expression *> props ;
		SparsePetriNet * spec;

		// form is a pair <isPlace, index>
		// isPlace false = it's a transition
		typedef std::pair<bool,int> node_t;

		typedef std::pair< std::pair<std::string, std::string>, int> arc_t;
		typedef std::vector<arc_t *> arcs_t;
		arcs_t topatch;

		std::string lastseen;
		bool dotext;
		bool isLTL;

		// map place name to index
		typedef ext_hash_map<std::string, int > index_t;
		index_t index;

	public :
		ExprHandler(SparsePetriNet * spec, bool isLTL=true):spec(spec),dotext(false),isLTL(isLTL) {

			for (int i=0, ie=spec->getPnames().size() ; i < ie ; ++i) {
				index_t::accessor acc;
				index.insert(acc, spec->getPnames().at(i));
				acc->second = i;
			}
		}

		static void characters(void * userData, const XML_Char * chars, int length) {
			ExprHandler * tthis = (ExprHandler *) userData;
			if (tthis->dotext) {
				tthis->lastseen = std::string(chars,length);
			}
		}


		/** {@inheritDoc} */
		static void startElement(void * userData, const XML_Char *name, const XML_Char **atts) {
			ExprHandler * tthis = (ExprHandler *) userData;
			std::string baliseName (name);
			// std::cerr << "start " << baliseName << std::endl;
			if ("property"== baliseName) { //$NON-NLS-1$
				Property * pdesc = new Property();
				tthis->stack.push(pdesc);
			} else if ("formula"== baliseName) { //$NON-NLS-1$
				// NOTHING
			} else if ("tokens-count"== baliseName) { //$NON-NLS-1$
				tthis->stack.push(new NaryOp(CARD));
			} else if ("place-bound"== baliseName) { //$NON-NLS-1$
				tthis->stack.push(new NaryOp(BOUND));
			} else if ("is-fireable"== baliseName) { //$NON-NLS-1$
				tthis->stack.push(new NaryOp(ENABLED));
			} else if ("property-set"== baliseName) { //$NON-NLS-1$
				// NOTHING
			} else if ("integer-le"== baliseName) { //$NON-NLS-1$
				// NOTHING
			} else if ("all-paths"== baliseName) { //$NON-NLS-1$
				// NOTHING
			} else if ("finally"== baliseName) { //$NON-NLS-1$
				// NOTHING
			} else if ("exists-path"== baliseName) { //$NON-NLS-1$
				// NOTHING
			} else if ("negation"== baliseName) { //$NON-NLS-1$
				// NOTHING
			} else if ("disjunction"== baliseName) { //$NON-NLS-1$
				// prepare Nary operator
				tthis->stack.push(nullptr);
			} else if ("conjunction"== baliseName) { //$NON-NLS-1$
				// prepare Nary operator
				tthis->stack.push(nullptr);
			} else if ("until"== baliseName) { //$NON-NLS-1$
				// NOTHING
			} else if ("before"== baliseName) { //$NON-NLS-1$
				// NOTHING
			} else if ("reach"== baliseName) { //$NON-NLS-1$
				// NOTHING
			} else if ("deadlock"== baliseName) {
				// NOTHING
			} else if ("next"== baliseName) { //$NON-NLS-1$
				// NOTHING
			} else if ("description"== baliseName) { //$NON-NLS-1$
				tthis->dotext = true;
			} else if ("place"== baliseName) { //$NON-NLS-1$
				tthis->dotext = true;
			} else if ("integer-constant"== baliseName) { //$NON-NLS-1$
				tthis->dotext = true;
			} else if ("id"== baliseName) { //$NON-NLS-1$
				tthis->dotext = true;
			} else if ("transition"== baliseName) { //$NON-NLS-1$
				tthis->dotext = true;
			} else if ("globally"== baliseName) { //$NON-NLS-1$
				// NOTHING
			} else {

				std::cerr << "Unknown XML tag in source file: " << baliseName; //$NON-NLS-1$
			}
		}


		/** {@inheritDoc} */
		static void endElement(void *userData, const XML_Char *name) {
			ExprHandler * tthis = (ExprHandler *) userData;
			std::string baliseName (name);
			// std::cerr << "end " << baliseName << std::endl;

			if ("property"== baliseName) { //$NON-NLS-1$
				tthis->spec->getProperties().push_back((Property *) tthis->stack.top());
				tthis->stack.pop();
			} else if ("formula"== baliseName) { //$NON-NLS-1$
				Expression * child = (Expression *) tthis->stack.top();
				tthis->stack.pop();
				Property * pdesc = (Property *) tthis->stack.top();
				pdesc->setBody(child);
			} else if ("integer-le"== baliseName) { //$NON-NLS-1$
				tthis->popBinary(LEQ);
			} else if ("negation"== baliseName) { //$NON-NLS-1$
				auto e = (Expression*) tthis->stack.top();
				tthis->stack.pop();
				tthis->stack.push(new BinOp(NOT,e,nullptr));
			} else if ("is-fireable"== baliseName) { //$NON-NLS-1$

			} else if ("tokens-count"== baliseName) { //$NON-NLS-1$

			} else if ("deadlock"== baliseName) {
				tthis->stack.push(new BinOp(DEAD, nullptr, nullptr));
			} else if ("description"== baliseName) { //$NON-NLS-1$
				std::string name = tthis->lastseen;
				// Property prop = (Property) stack.peek();
				// prop.setComment(name);
				tthis->dotext = false;

			} else if ("place"== baliseName) {
				NaryOp * context = (NaryOp *) tthis->stack.top();
				context->addChild(new VarRef(tthis->findPlace(tthis->lastseen)));
				tthis->dotext = false;

			} else if ("integer-constant"== baliseName) {
				auto lastint = std::stol(tthis->lastseen);
				tthis->stack.push(new Constant(lastint));
				tthis->dotext = false;

			} else if ("id"== baliseName) { //$NON-NLS-1$
				Property *prop = (Property*) tthis->stack.top();
				prop->setName(tthis->lastseen);
				tthis->dotext = false;

			} else if ("disjunction"== baliseName) { //$NON-NLS-1$
				tthis->popNary(OR);
			} else if ("conjunction"== baliseName) { //$NON-NLS-1$
				tthis->popNary(AND);
			} else if ("transition"== baliseName) {
				//NaryOp* enab = (NaryOp*) tthis->stack.top();
				//enab->addChild(Expression.trans(findTransition(name)));
				throw "Enabling criterion not supported currently.";
				tthis->dotext = false;
			} else if ("before" == baliseName) { //$NON-NLS-1$
				// NOTHING
			} else if ("reach" == baliseName) { //$NON-NLS-1$
				// NOTHING
			} else if (! tthis->isLTL) {
				// temporal operator handling for CTL properties
				if ( ("globally"== baliseName || "finally"== baliseName || "next"== baliseName || "until"== baliseName ) ) {
					tthis->stack.push(new std::string(baliseName));
				} else if ("all-paths"== baliseName) { //$NON-NLS-1$
					std::string * childbalise = (std::string *) tthis->stack.top();
					tthis->stack.pop();
					if (*childbalise=="globally") {
						auto e = (Expression*) tthis->stack.top();
						tthis->stack.pop();
						tthis->stack.push(new BinOp(AG, e, nullptr));
					} else if (*childbalise=="finally") {
						auto e = (Expression*) tthis->stack.top();
						tthis->stack.pop();tthis->stack.push(new BinOp(AF, e, nullptr));
					} else if (*childbalise=="next") {
						auto e = (Expression*) tthis->stack.top();
						tthis->stack.pop();
						tthis->stack.push(new BinOp(AX, e, nullptr));
					} else if (*childbalise=="until") {
						tthis->popBinary(AU);
					}
				} else if ("exists-path"== baliseName) { //$NON-NLS-1$
					std::string * childbalise = (std::string *) tthis->stack.top();
					tthis->stack.pop();
					if (*childbalise=="globally") {
						auto e = (Expression*) tthis->stack.top();
						tthis->stack.pop();
						tthis->stack.push(new BinOp(EG, e, nullptr));
					} else if (*childbalise=="finally") {
						auto e = (Expression*) tthis->stack.top();
						tthis->stack.pop();tthis->stack.push(new BinOp(EF, e, nullptr));
					} else if (*childbalise=="next") {
						auto e = (Expression*) tthis->stack.top();
						tthis->stack.pop();
						tthis->stack.push(new BinOp(EX, e, nullptr));
					} else if (*childbalise=="until") {
						tthis->popBinary(EU);
					}
				} else {
					std::cerr << "Unknown XML tag in source file: " << baliseName << std::endl; //$NON-NLS-1$
				}
			} else {
				// temporal operator handling for LTL properties
				if ("all-paths"== baliseName) {
					// only as first node
					// hence leave stack alone and skip this node
				} else if ("globally"== baliseName) {
					auto e = (Expression*) tthis->stack.top();
					tthis->stack.pop();
					tthis->stack.push(new BinOp(G, e, nullptr));
				} else if ("finally"== baliseName) {
					auto e = (Expression*) tthis->stack.top();
					tthis->stack.pop();
					tthis->stack.push(new BinOp(F, e, nullptr));
				} else if ("next"== baliseName) {
					auto e = (Expression*) tthis->stack.top();
					tthis->stack.pop();
					tthis->stack.push(new BinOp(X, e, nullptr));
				} else if ("until"== baliseName) {
					tthis->popBinary(U);
				} else {
					std::cerr << "Unknown XML tag in source file: " << baliseName << std::endl; //$NON-NLS-1$
				}
			}
		}


		void popNary(Op op) {
			std::vector<Expression* > operands;
			while (stack.top() != nullptr) {
				Expression * r = (Expression*) stack.top();
				stack.pop();
				operands.push_back(r);
			}
			// the nullptr separator
			stack.pop();
			stack.push(new NaryOp(op, operands));
		}


		void popBinary(Op op) {
			Expression * r = (Expression *) stack.top();
			stack.pop();
			Expression * l = (Expression *) stack.top();
			stack.pop();
			stack.push(new BinOp(op, l, r));
		}

		int findPlace (const std::string & name) {
			index_t::const_accessor acc;
			index.find(acc, name);
			return acc->second;
		}

		/**
		 * @return the order loaded from the XML file
		 */
		std::vector<Expression *> getParseResult() {
			return props;
		}
	};

}

#endif
