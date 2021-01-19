/*
 * PTNetLoader.cpp
 *
 * Created on: Oct 13, 2020
 * Author: ythierry
 */

#include "expr/parse/ExprLoader.h"

#include <cstdio>
#include <expat.h>
#include "expr/parse/ExprHandler.h"

namespace petri::expr {

	int loadXMLProperties(std::string filename, SparsePetriNet * spec) {
		char buf[BUFSIZ];
		XML_Parser parser = XML_ParserCreate(NULL);
		int done;
		ExprHandler handler(spec,true);

		FILE * in = fopen(filename.c_str(),"r");

		auto initial = spec->getProperties().size();
		XML_SetUserData(parser, &handler);
		XML_SetElementHandler(parser, &ExprHandler::startElement, &ExprHandler::endElement);
		XML_SetCharacterDataHandler(parser, &ExprHandler::characters);
		do {
			size_t len = fread(buf, 1, sizeof(buf), in);
			done = len < sizeof(buf);
			if (XML_Parse(parser, buf, (int)len, done) == XML_STATUS_ERROR)
			{
				fprintf(stderr, "%s at line %d \n",
						XML_ErrorString(XML_GetErrorCode(parser)),
						(int)XML_GetCurrentLineNumber(parser));
				XML_ParserFree(parser);
				return -1;
			}
		} while (! done);
		XML_ParserFree(parser);
		fclose(in);
		return spec->getProperties().size() - initial;
	}

}
