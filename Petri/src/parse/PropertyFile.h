/*
 * PropertyFile.h
 *
 * Choice of the property syntax: MCC XML or s-expressions, by explicit
 * request or by file extension (".xml" is MCC, anything else s-expressions).
 */
#ifndef PETRI_PARSE_PROPERTYFILE_H_
#define PETRI_PARSE_PROPERTYFILE_H_

#include <string>
#include <vector>

#include "core/SparsePetriNet.h"
#include "expr/Property.h"
#include "parse/mcc/PropertyLoader.h"
#include "parse/sexpr/PropertyReader.h"

namespace petri
{

enum class PropertySyntax
{
  Auto, Mcc, Sexpr
};

/** Parse "auto" | "mcc" | "sexpr"; throws std::string otherwise. */
inline PropertySyntax propertySyntaxOf (const std::string &name)
{
  if (name == "auto") return PropertySyntax::Auto;
  if (name == "mcc") return PropertySyntax::Mcc;
  if (name == "sexpr") return PropertySyntax::Sexpr;
  throw std::string ("Unknown property syntax: " + name + " (auto, mcc, sexpr)");
}

inline PropertySyntax syntaxByExtension (const std::string &filename)
{
  size_t dot = filename.rfind ('.');
  std::string ext = dot == std::string::npos ? "" : filename.substr (dot);
  return ext == ".xml" || ext == ".XML" ? PropertySyntax::Mcc : PropertySyntax::Sexpr;
}

/** Load a property file in the given or detected syntax. Throws std::string. */
template<typename T>
  std::vector<petri::expr::Property> loadPropertyFile (const std::string &filename,
                                                       const SparsePetriNet<T> &net,
                                                       PropertySyntax syntax = PropertySyntax::Auto)
  {
    if (syntax == PropertySyntax::Auto) syntax = syntaxByExtension (filename);
    if (syntax == PropertySyntax::Mcc) return petri::mcc::loadProperties (filename, net);
    return petri::sexpr::loadProperties (filename, net);
  }

} // namespace petri

#endif /* PETRI_PARSE_PROPERTYFILE_H_ */
