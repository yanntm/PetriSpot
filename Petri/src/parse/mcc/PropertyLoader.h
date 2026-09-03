/*
 * PropertyLoader.h
 *
 * Drive expat over an MCC property XML file and return the properties it
 * declares, resolved against a net.
 */
#ifndef PETRI_PARSE_MCC_PROPERTYLOADER_H_
#define PETRI_PARSE_MCC_PROPERTYLOADER_H_

#include <chrono>
#include <cstdio>
#include <expat.h>
#include <string>
#include <vector>

#include "core/Log.h"
#include "core/SparsePetriNet.h"
#include "expr/Property.h"
#include "parse/mcc/PropertyHandler.h"

namespace petri::mcc
{

/** Parse filename; throws a std::string on an unresolved place/transition. */
template<typename T>
  std::vector<petri::expr::Property> loadProperties (const std::string &filename,
                                                     const SparsePetriNet<T> &net)
  {
    auto time = std::chrono::steady_clock::now ();
    petri::writeToLog ("Parsing property file : " + filename);

    FILE *in = fopen (filename.c_str (), "r");
    if (!in) {
      throw std::string ("Cannot open property file: " + filename);
    }
    PropertyHandler<T> handler (net);
    XML_Parser parser = XML_ParserCreate (NULL);
    XML_SetUserData (parser, &handler);
    XML_SetElementHandler (parser, &PropertyHandler<T>::startElement,
                           &PropertyHandler<T>::endElement);
    XML_SetCharacterDataHandler (parser, &PropertyHandler<T>::characters);

    char buf[BUFSIZ];
    bool done = false;
    while (!done) {
      size_t len = fread (buf, 1, sizeof(buf), in);
      done = len < sizeof(buf);
      if (XML_Parse (parser, buf, static_cast<int> (len), done) == XML_STATUS_ERROR) {
        std::string err = std::string (XML_ErrorString (XML_GetErrorCode (parser)))
            + " at line " + std::to_string (XML_GetCurrentLineNumber (parser));
        XML_ParserFree (parser);
        fclose (in);
        throw err;
      }
    }
    XML_ParserFree (parser);
    fclose (in);

    std::vector<petri::expr::Property> result = std::move (handler.getParseResult ());
    petri::writeToLog (
        "Parsed " + std::to_string (result.size ()) + " properties in "
            + std::to_string (
                std::chrono::duration_cast<std::chrono::milliseconds> (
                    std::chrono::steady_clock::now () - time).count ()) + " ms.");
    return result;
  }

} // namespace petri::mcc

#endif /* PETRI_PARSE_MCC_PROPERTYLOADER_H_ */
