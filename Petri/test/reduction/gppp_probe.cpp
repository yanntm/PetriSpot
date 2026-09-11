#include <iostream>
#include <memory>
#include "parse/PTNetLoader.h"
#include "parse/PropertyFile.h"
#include "cli/CtlDriver.h"

int main(int argc, char **argv) {
  if (argc != 2) return 2;
  std::string dir = argv[1];
  std::unique_ptr<SparsePetriNet<long>> net(loadXML<long>(dir + "/model.pnml"));
  auto props = petri::loadPropertyFile(dir + "/CTLCardinality.xml", *net);
  auto f = petri::expr::ctlNormalize(props.at(8).ctl);
  petri::walk::WalkNet<long> wn(*net);
  petri::cli::Options options;
  options.ctlRounds = 1; options.ctlSteps = 100; options.ctlRegion = 100;
  options.ctlRunLength = 100;
  auto outcome = petri::cli::checkCtl(options, wn, f.kids.at(1), 1 + 7919 * 8, 1000);
  std::cout << "AG verdict " << petri::ctl::to_string(outcome.verdict) << '\n';
  outcome.evidence.print(std::cout, &net->getTnames());
  petri::ctl::Cursor<long> cur(wn);
  for (auto t : outcome.evidence.path) {
    if (!cur.enabled.isEnabled(t)) throw std::runtime_error("Disabled evidence transition");
    cur.fire(t);
  }
  for (size_t p = 0; p < net->getPlaceCount(); ++p)
    std::cout << net->getPnames()[p] << '=' << cur.marking.get(p) << '\n';
}
