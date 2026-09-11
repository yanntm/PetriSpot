#pragma once
#include <filesystem>
#include <fstream>
#include <memory>
#include "cli11/CLI11.hpp"
#include "expr/SexprPrinter.h"
#include "io/PNETIO.h"
#include "parse/PTNetLoader.h"
#include "parse/PropertyFile.h"
#include "reduction/Pipeline.h"

namespace petri::reduction {
/** Standalone transformation: publish a paired PNET/property directory, without
 * invoking a solver. PNET indices and emitted formula indices agree; names are
 * exported separately because PNET's mandatory payload has no name table. */
template<class T> int reduceCommand(int argc, char** argv) {
  std::string pnml, pnet, propertyFile, syntax = "auto", output;
  long query = -1, milliseconds = 15000;
  bool noAgglo = false;
  CLI::App app("Reduce a P/T net and its formulas; write a matched pair without solving.");
  auto* input = app.add_option("-i", pnml, "Input PNML net.");
  app.add_option("--net", pnet, "Input PNET net.")->excludes(input);
  app.add_option("--props", propertyFile, "MCC XML or s-expression formulas.")->required();
  app.add_option("--propsSyntax", syntax, "auto, mcc or sexpr.");
  app.add_option("--query", query, "Select one formula by zero-based index.");
  app.add_option("--output", output, "New directory for model.pnet, properties.sexpr and names.sexpr.")->required();
  app.add_option("--reductionMs", milliseconds, "Reduction time limit in milliseconds.")->check(CLI::PositiveNumber);
  app.add_flag("--reductionNoAgglo", noAgglo, "Disable agglomeration.");
  try { app.parse(argc, argv); } catch (const CLI::ParseError& error) { return app.exit(error); }
  if (pnml.empty() && pnet.empty()) throw std::invalid_argument("Specify -i or --net");
  namespace fs = std::filesystem;
  fs::path destination(output);
  if (destination.filename().empty()) destination = destination.parent_path();
  fs::path staging = destination.string() + ".tmp";
  if (fs::exists(destination) || fs::exists(staging))
    throw std::invalid_argument("Reduction output or staging directory already exists");

  typename PNETIO<T>::Blocks blocks;
  std::unique_ptr<SparsePetriNet<T>> original(pnet.empty() ? loadXML<T>(pnml) : PNETIO<T>::read(pnet, &blocks));
  if (!original) throw std::runtime_error("Cannot read reduction input");
  if (!blocks.empty()) std::cerr << "Reduction export omits input counting records; their transport is not implemented.\n";
  auto properties = petri::loadPropertyFile<T>(propertyFile, *original, petri::propertySyntaxOf(syntax));
  if (query < -1 || (query >= 0 && static_cast<size_t>(query) >= properties.size()))
    throw std::invalid_argument("Reduction query index is out of range");
  if (query >= 0) {
    auto selected = properties[static_cast<size_t>(query)]; properties.clear(); properties.push_back(std::move(selected));
  }
  if (properties.empty()) throw std::invalid_argument("No formulas to reduce");
  for (const auto& property : properties)
    if (property.kind == expr::PropertyKind::Unsupported || !property.hint.empty())
      throw std::invalid_argument("Reduction export needs supported formulas without transition hints");
  Configuration config;
  config.timeLimit = std::chrono::milliseconds(milliseconds); config.agglomeration = !noAgglo;
  auto prepared = prepare(*original, properties, true, false, config, nullptr, std::cerr);
  const auto* reduced = &prepared;

  if (!fs::create_directory(staging)) throw std::runtime_error("Cannot create reduction staging directory");
  try {
    std::ofstream net(staging / "model.pnet", std::ios::binary);
    if (!PNETIO<T>::write(reduced->net, net)) throw std::runtime_error("Cannot write reduced PNET");
    net.close();
    if (!net) throw std::runtime_error("Cannot close reduced PNET");
    std::ofstream formulas(staging / "properties.sexpr");
    for (const auto& property : properties) { expr::printSexpr(formulas, property, nullptr); formulas << '\n'; }
    formulas.close();
    if (!formulas) throw std::runtime_error("Cannot write reduced formulas");
    std::ofstream names(staging / "names.sexpr");
    for (size_t p = 0; p < reduced->net.getPlaceCount(); ++p) {
      names << "(place p" << p << " "; expr::printSexprName(names, reduced->net.getPnames()[p]); names << ")\n";
    }
    for (size_t t = 0; t < reduced->net.getTransitionCount(); ++t) {
      names << "(transition t" << t << " "; expr::printSexprName(names, reduced->net.getTnames()[t]); names << ")\n";
    }
    names.close();
    if (!names) throw std::runtime_error("Cannot write reduction names");
    fs::rename(staging, destination);
  } catch (...) {
    std::error_code ignored; fs::remove_all(staging, ignored); throw;
  }
  std::cout << "Reduced model: " << (destination / "model.pnet").string()
      << "\nReduced formulas: " << (destination / "properties.sexpr").string() << '\n';
  return 0;
}
}
