//
// Filter events based on a good run list
//

#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/DbService/inc/GrlList.hh"
#include "Offline/DbService/inc/GrlTool.hh"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include <iostream>
#include <memory>
#include <vector>

using namespace std;

namespace mu2e {

class GoodRunFilter : public art::EDFilter {
public:
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::OptionalAtom<string> name{Name("name"), Comment("name of list in database")};
    fhicl::OptionalAtom<string> fileName{Name("fileName"), Comment("name of file containing a list")};
  };
  // typedef art::EDFilter::Table<Config> Parameters;
  using Parameters = art::EDFilter::Table<Config>;
  explicit GoodRunFilter(const Parameters& conf);
  virtual ~GoodRunFilter() {}

  bool filter(art::Event& event);

private:
  std::unique_ptr<GrlList> _list;
};

GoodRunFilter::GoodRunFilter(const Parameters& conf) : art::EDFilter{conf} {

  string name;
  if (conf().name()) {
    conf().name(name);
    GrlTool tool;
    _list = make_unique<GrlList>(tool.list(name));
  } else if (conf().fileName()) {
    conf().fileName(name);
    ConfigFileLookupPolicy configFile;
    std::string fn = configFile(name);
    GrlHeader gh("fromFile");
    _list = make_unique<GrlList>(gh, fn);
  } else {
    throw cet::exception("GRL_FILTER_CONFIG") << "GoodRunFilter not configured\n";
  }
}

bool GoodRunFilter::filter(art::Event& event) {
  return _list->goodSubRun(event.run(), event.subRun());
}

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::GoodRunFilter)
