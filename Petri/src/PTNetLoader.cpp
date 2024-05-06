/*
 * PTNetLoader.cpp
 *
 *  Created on: Oct 13, 2020
 *      Author: ythierry
 */

#include "PTNetLoader.h"
#include "PTNetHandler.h"
#include "InvariantMiddle.h"
#include <cstdio>
#include <expat.h>


SparsePetriNet * loadXML(std::string filename) {
	char buf[BUFSIZ];
	XML_Parser parser = XML_ParserCreate(NULL);
	int done;
	PTNetHandler handler;

	FILE * in = fopen(filename.c_str(),"r");

	std::string logMessage = "Parsing pnml file : ";
	logMessage += filename;

	InvariantMiddle::writeToLog(logMessage);

	XML_SetUserData(parser, &handler);
	XML_SetElementHandler(parser, &PTNetHandler::startElement, &PTNetHandler::endElement);
	XML_SetCharacterDataHandler(parser, &PTNetHandler::characters);
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
	return handler.getParseResult();
}


