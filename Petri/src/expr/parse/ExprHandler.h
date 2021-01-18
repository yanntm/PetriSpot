#ifndef PTNETHANDLER_H_
#define PTNETHANDLER_H_

#include "SparsePetriNet.h"
#include "expr/NaryOp.h"
#include "expr/BinOp.h"
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
		const SparsePetriNet * spec;

		// form is a pair <isPlace, index>
		// isPlace false = it's a transition
		typedef std::pair<bool,int> node_t;

		// map object name to <isPlace,index>
		typedef ext_hash_map<std::string, node_t > index_t;
		index_t index;

		typedef std::pair< std::pair<std::string, std::string>, int> arc_t;
		typedef std::vector<arc_t *> arcs_t;
		arcs_t topatch;

		std::string lastseen;
		bool dotext;
		bool isLTL;

	public :
		ExprHandler(const SparsePetriNet * spec, bool isLTL=true):spec(spec),dotext(false),isLTL(isLTL) {}

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
			if ("property"== baliseName) { //$NON-NLS-1$
				//Property pdesc = new Property();
				// stack.push(pdesc);
			} else if ("formula"== baliseName) { //$NON-NLS-1$
				// NOTHING
			} else if ("tokens-count"== baliseName) { //$NON-NLS-1$
				stack.push(new NaryOp(CARD));
			} else if ("place-bound"== baliseName) { //$NON-NLS-1$
				stack.push(new NaryOp(BOUND));
			} else if ("is-fireable"== baliseName) { //$NON-NLS-1$
				stack.push(new NaryOp(ENABLED));
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
				stack.push(nullptr);
			} else if ("conjunction"== baliseName) { //$NON-NLS-1$
				// prepare Nary operator
				stack.push(nullptr);
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
				dotext = true;
			} else if ("place"== baliseName) { //$NON-NLS-1$
				dotext = true;
			} else if ("integer-constant"== baliseName) { //$NON-NLS-1$
				dotext = true;
			} else if ("id"== baliseName) { //$NON-NLS-1$
				dotext = true;
			} else if ("transition"== baliseName) { //$NON-NLS-1$
				dotext = true;
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
			if ("property"== baliseName) { //$NON-NLS-1$
				// spec.getProperties().add((Property) stack.pop());
			} else if ("formula"== baliseName) { //$NON-NLS-1$
				Expression * child = (Expression *) stack.pop();
				// Property pdesc = (Property) stack.peek();
				// pdesc.setBody(child);
				props.push_back(child);

			} else if ("integer-le"== baliseName) { //$NON-NLS-1$
				tthis->popBinary(LEQ);
			} else if ("negation"== baliseName) { //$NON-NLS-1$
				stack.push(new BinOp(NOT,(Expression*) stack.pop(),nullptr));
			} else if ("is-fireable"== baliseName) { //$NON-NLS-1$

			} else if ("tokens-count"== baliseName) { //$NON-NLS-1$

			} else if ("deadlock"== baliseName) {
				stack.push(new BinOp(DEAD, nullptr, nullptr));
			} else if ("description"== baliseName) { //$NON-NLS-1$
				std::string name = lastseen;
				// Property prop = (Property) stack.peek();
				// prop.setComment(name);
				dotext = false;

			} else if ("place"== baliseName) {
				String name = (String) stack.pop();
				NaryOp * context = (NaryOp *) stack.top();
				context->addChild(new VarRef(findPlace(name)));
				dotext = false;

			} else if ("integer-constant"== baliseName) {
				String name = (String) stack.pop();
				stack.push(Expression.constant(Integer.parseInt(name)));
				dotext = false;

			} else if ("id"== baliseName) { //$NON-NLS-1$
				String name = (String) stack.pop();
				Property prop = (Property) stack.peek();
				prop.setName(name);
				dotext = false;

			} else if ("disjunction"== baliseName) { //$NON-NLS-1$
				popNary(Op.OR);
			} else if ("conjunction"== baliseName) { //$NON-NLS-1$
				popNary(Op.AND);
			} else if ("transition"== baliseName) {
				String name = (String) stack.pop();
				NaryOp enab = (NaryOp) stack.peek();
				enab.addChild(Expression.trans(findTransition(name)));
				dotext = false;
			} else if (! isLTL) {
				// temporal operator handling for CTL properties
				if ( ("globally"== baliseName || "finally"== baliseName || "next"== baliseName || "until"== baliseName ) ) {
					stack.push(baliseName);
				} else if ("all-paths"== baliseName) { //$NON-NLS-1$
					String childbalise = (String) stack.pop();
					if (childbalise.equals("globally")) {
						stack.push(Expression.op(Op.AG, (Expression) stack.pop(), null));
					} else if (childbalise.equals("finally")) {
						stack.push(Expression.op(Op.AF, (Expression) stack.pop(), null));
					} else if (childbalise.equals("next")) {
						stack.push(Expression.op(Op.AX, (Expression) stack.pop(), null));
					} else if (childbalise.equals("until")) {
						popBinary(Op.AU);
					}
				} else if ("exists-path"== baliseName) { //$NON-NLS-1$
					String childbalise = (String) stack.pop();
					if (childbalise.equals("globally")) {
						stack.push(Expression.op(Op.EG, (Expression) stack.pop(), null));
					} else if (childbalise.equals("finally")) {
						stack.push(Expression.op(Op.EF, (Expression) stack.pop(), null));
					} else if (childbalise.equals("next")) {
						stack.push(Expression.op(Op.EX, (Expression) stack.pop(), null));
					} else if (childbalise.equals("until")) {
						popBinary(Op.EU);
					}
				}
			} else {
				// temporal operator handling for LTL properties
				if ("all-paths"== baliseName) {
					// only as first node
					// hence leave stack alone and skip this node
				} else if ("globally"== baliseName) {
					stack.push(Expression.op(Op.G, (Expression) stack.pop(), null));
				} else if ("finally"== baliseName) {
					stack.push(Expression.op(Op.F, (Expression) stack.pop(), null));
				} else if ("next"== baliseName) {
					stack.push(Expression.op(Op.X, (Expression) stack.pop(), null));
				} else if ("until"== baliseName) {
					popBinary(Op.U);
				}
			}


			// Balise MODEL
			if ("toolspecific"==baliseName) {
				tthis->inOpaqueToolSpecific = false;
			} else if (tthis->inOpaqueToolSpecific) {
				// skipping this stuff
				return;
			} else if ("net"==baliseName) { //$NON-NLS-1$
				tthis->stack.pop();
				assert(tthis->stack.empty());
			} else if ("name"==baliseName) { //$NON-NLS-1$
				// names of objects are dropped, we only use identifiers.
				tthis->readtext = false;
				tthis->lastseen = "";
			} else if ("page"==baliseName) { //$NON-NLS-1$
				// ignored pages
			} else if ("place"==baliseName) {
				tthis->stack.pop();
			} else if ("transition"==baliseName) {
				tthis->stack.pop();
			} else if ("arc"==baliseName) {
				arc_t * arc = (arc_t*)tthis->stack.top();

				index_t::const_accessor accsrc;
				bool oksrc = tthis->index.find(accsrc,arc->first.first);

				index_t::const_accessor acctgt;
				bool oktgt = tthis->index.find(acctgt,arc->first.second);

				if (oktgt && oksrc) {
					if (accsrc->second.first) {
						// source is a place
						tthis->net->addPreArc(accsrc->second.second, acctgt->second.second, arc->second);
					} else {
						tthis->net->addPostArc(acctgt->second.second, accsrc->second.second, arc->second);
					}
					delete arc;
				} else {
					tthis->topatch.push_back(arc);
				}
				tthis->stack.pop();
			} else if ("text"==baliseName) {
				tthis->doIt = false;
			} else if ("initialMarking"==baliseName) {
				size_t p = (size_t) tthis->stack.top();
				tthis->net->setMarking(p,tthis->lastint);
				tthis->readint = false;
				tthis->lastint = -1;
			} else if ("inscription"==baliseName) {
				arc_t * arc = (arc_t *) tthis->stack.top();

				arc->second = tthis->lastint;
				tthis->readint = false;
				tthis->lastint = -1;
			} else if ("graphics"==baliseName || "offset"==baliseName || "position"==baliseName || "fill"==baliseName || "line"==baliseName || "dimension"==baliseName) {
				//skip
			} else if ("pnml"==baliseName) {
				// patch missing arc targets
				for (arc_t * & arc : tthis->topatch) {
					index_t::const_accessor accsrc;
					bool oksrc = tthis->index.find(accsrc,arc->first.first);

					index_t::const_accessor acctgt;
					bool oktgt = tthis->index.find(acctgt,arc->first.second);

					if (oktgt && oksrc) {
						if (accsrc->second.first) {
							// source is a place
							tthis->net->addPreArc(accsrc->second.second, acctgt->second.second, arc->second);
						} else {
							tthis->net->addPostArc(acctgt->second.second, accsrc->second.second, arc->second);
						}
						delete arc;
					} else {
						std::string  err = "Problem when linking arc : source or target node not found <" + accsrc->first + "," + acctgt->first + ">";
						throw err.c_str();
					}
					delete arc;
				}
				tthis->topatch.clear();
			} else {
				std::cerr << "Unknown XML tag in source file: " << baliseName << std::endl; //$NON-NLS-1$
			}
		}


		void popNary(Op op) {
			std::vector<Expression* > operands;
			while (stack.top() != nullptr) {
				Expression * r = (Expression*) stack.pop();
				operands.push_back(r);
			}
			// the nullptr separator
			stack.pop();
			stack.push(new NaryOp(op, operands));
		}


		void popBinary(Op op) {
			Expression * r = (Expression *) stack.pop();
			Expression * l = (Expression *) stack.pop();
			stack.push(new BinOp(op, l, r));
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
