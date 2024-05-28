/*
 * PTNetLoader.h
 *
 *  Created on: Oct 13, 2020
 *      Author: ythierry
 */

#ifndef PTNETLOADER_H_
#define PTNETLOADER_H_

#include "SparsePetriNet.h"
#include "PTNetHandler.h"
#include "InvariantMiddle.h"
#include <cstdio>
#include <expat.h>


template<typename T>
SparsePetriNet<T> * loadXML(std::string filename) {
	char buf[BUFSIZ];
	XML_Parser parser = XML_ParserCreate(NULL);
	int done;
	PTNetHandler<T> handler;
	auto time = std::chrono::steady_clock::now();

	FILE * in = fopen(filename.c_str(),"r");

	std::string logMessage = "Parsing pnml file : " + filename;
	InvariantMiddle<T>::writeToLog(logMessage);

	XML_SetUserData(parser, &handler);
	XML_SetElementHandler(parser, &PTNetHandler<T>::startElement, &PTNetHandler<T>::endElement);
	XML_SetCharacterDataHandler(parser, &PTNetHandler<T>::characters);
	do {
		size_t len = fread(buf, 1, sizeof(buf), in);
		done = len < sizeof(buf);
		if (XML_Parse(parser, buf, (int)len, done) == XML_STATUS_ERROR) {
			fprintf(stderr, "%s at line %ld \n",
					XML_ErrorString(XML_GetErrorCode(parser)),
					XML_GetCurrentLineNumber(parser));
			XML_ParserFree(parser);
			return nullptr;
		}
	} while (! done);
	XML_ParserFree(parser);
	fclose(in);

	logMessage = "Parsed PT model containing " + std::to_string(handler.getParseResult()->getPlaceCount())
			+ " places and " + std::to_string(handler.getParseResult()->getTransitionCount())
			+ " transitions and " + std::to_string(handler.getParseResult()->getArcCount())
			+ " arcs in " + std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>
			(std::chrono::steady_clock::now() - time).count()) + " ms.";
	InvariantMiddle<T>::writeToLog(logMessage);

	return handler.getParseResult();
}


#endif /* PTNETLOADER_H_ */
