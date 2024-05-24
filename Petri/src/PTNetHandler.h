#ifndef PTNETHANDLER_H_
#define PTNETHANDLER_H_

#include "SparsePetriNet.h"
#include "ddd/util/ext_hash_map.hh"

#include <expat.h>
#include <stack>

/**
 * A class to parse a PT model from an PNML file.
 * @author Yann Thierry-Mieg 2015
 */
class PTNetHandler {

	// context stack
	std::stack<void *> stack;

	// object constructed
	SparsePetriNet * net = new SparsePetriNet();

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
	bool readtext;

	long lastint;
	bool readint;
	bool inOpaqueToolSpecific;

	bool doIt;

public :
	 PTNetHandler():net(new SparsePetriNet()),readtext(false),lastint(-1),readint(false),inOpaqueToolSpecific(false),doIt(false) {}

	 static void characters(void * userData, const XML_Char * chars, int length) {
		 PTNetHandler * tthis = (PTNetHandler *) userData;
		 if (tthis->inOpaqueToolSpecific) {
			return;
		} else if (tthis->doIt) {
			if (tthis->readtext) {
				tthis->lastseen = std::string(chars,length);
			} else if (tthis->readint) {
				std::string laststr = std::string(chars,length);
				tthis->lastint = std::stol(laststr);
			}
		} 
	}
	
	
	/** {@inheritDoc} */
	 static void startElement(void * userData, const XML_Char *name, const XML_Char **atts) {
//		if (doNupn) {
//			nupnHandler.startElement(uri, localName, baliseName, attributes);
//		} else
		PTNetHandler * tthis = (PTNetHandler *) userData;
		std::string baliseName (name);
		if ("net" == baliseName) { //$NON-NLS-1$

			for (int i=0 ; atts[i] != nullptr ; i+=2) {
				std::string bname = atts[i];
				std::string bval = atts[i+1];

				if ("id" == bname) {
					tthis->net->setName(bval);
				} else if ("type" == bname) {
					if ("http://www.pnml.org/version-2009/grammar/ptnet" != bval) {
						throw "Net is not a P/T net-> Colors are not supported currently.";
					}
				}
			}

			tthis->stack.push(& tthis->net);
		} else if ("name"==baliseName) {
			tthis->readtext = true;
			
		} else if ("page"==baliseName) {
			//SparsePetriNet * pn = (SparsePetriNet *) stack.top();
			// pages are ignored currently
//			Page page = PtnetFactory.eINSTANCE.createPage();
//			page.setId(attributes.getValue("id"));
//			pn.getPages().add(page);
//			stack.push(page);
		} else if ("place"==baliseName) {
			SparsePetriNet * pn =  tthis->net;

			std::string id;
			for (int i=0 ; atts[i] != nullptr ; i+=2) {
				std::string bname = atts[i];
				std::string bval = atts[i+1];

				if ("id" == bname) {
					id = bval;
					break;
				}
			}
			size_t pid = pn->addPlace(id, 0);
			index_t::accessor acc;
			tthis->index.insert(acc, id);
			acc->second = {true,pid};
			tthis->stack.push((void*)pid);
		} else if ("initialMarking"==baliseName) {
			tthis->readint = true;
		} else if ("inscription"==baliseName) {
			tthis->readint = true;
			
		} else if ("transition"==baliseName) {
			SparsePetriNet * pn = tthis->net;

			std::string id;
			for (int i=0 ; atts[i] != nullptr ; i+=2) {
				std::string bname = atts[i];
				std::string bval = atts[i+1];

				if ("id" == bname) {
					id = bval;
					break;
				}
			}
			size_t tid = pn->addTransition(id);
			index_t::accessor acc;
			tthis->index.insert(acc, id);
			acc->second = {false,tid};
			tthis->stack.push((void*)tid);
		} else if ("arc"==baliseName) {
			SparsePetriNet * pn = (SparsePetriNet *) tthis->stack.top();

			std::string source;
			std::string target;
			for (int i=0 ; atts[i] != nullptr ; i+=2) {
				std::string bname = atts[i];
				std::string bval = atts[i+1];

				if ("source" == bname) {
					source = bval;
				} else if ("target" == bname) {
					target = bval;
				}
				// ignore arc ID
			}

			// default arc weight is 1
			arc_t * arc = new arc_t({source,target},1);

			tthis->stack.push(arc);
		} else if ("toolspecific"==baliseName) {
			tthis->inOpaqueToolSpecific = true;
		} else if ("text"==baliseName) {
			tthis->doIt  = true;
		} else if ("graphics"==baliseName || "offset"==baliseName || "position"==baliseName || "fill"==baliseName || "line"==baliseName || "dimension"==baliseName) {
			//skip
		} else if ("pnml"==baliseName) {
			// skip
		} else {
			std::cerr << "Unknown XML tag in source file: " << baliseName; //$NON-NLS-1$
		}
	}


	/** {@inheritDoc} */
	 static void endElement(void *userData, const XML_Char *name) {
		 PTNetHandler * tthis = (PTNetHandler *) userData;
 		std::string baliseName (name);
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
					delete arc;
					throw err.c_str();
				}
			}
			tthis->topatch.clear();
		} else {
			std::cerr << "Unknown XML tag in source file: " << baliseName << std::endl; //$NON-NLS-1$
		}
	}




	/**
	 * @return the order loaded from the XML file
	 */
	 SparsePetriNet * getParseResult() {
		return net;
	}
};

#endif
