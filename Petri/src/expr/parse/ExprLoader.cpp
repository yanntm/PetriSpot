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

std::vector<Expression *> loadXML(std::string filename) {
  char buf[BUFSIZ];
  XML_Parser parser = XML_ParserCreate(NULL);
  int done;
  ExprHandler handler;

  FILE * in = fopen(filename.c_str(),"r");

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
	return nullptr;
      }
  } while (! done);
  XML_ParserFree(parser);
  fclose(in);
  return handler.getParseResult();
}

}
