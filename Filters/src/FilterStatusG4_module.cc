// Pass perfect events, fails events with StatusG4 > 0.  This allows
// to e.g. select events with trapped particles (by negating the
// filter) to study them.
//
// Andrei Gaponenko, 2013

#include <string>
#include <map>
#include <sstream>

// art includes.
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "MCDataProducts/inc/StatusG4.hh"

namespace mu2e {

  //================================================================
  class FilterStatusG4 : public art::EDFilter {
    art::InputTag input_;
    typedef std::map<int,int> StatMap;
    StatMap stats_;
  public:
    explicit FilterStatusG4(const fhicl::ParameterSet& pset);
    virtual bool filter(art::Event& event) override;
    virtual void endJob() override;
  };

  //================================================================
  FilterStatusG4::FilterStatusG4(const fhicl::ParameterSet& pset)
    : input_(pset.get<std::string>("input"))
  {}

  //================================================================
  bool FilterStatusG4::filter(art::Event& event) {
    auto ih = event.getValidHandle<StatusG4>(input_);
    ++stats_[ih->status()];
    return (ih->status() == 0);
  }

  //================================================================
  void FilterStatusG4::endJob() {
    std::ostringstream os;
    os<< "FilterStatusG4_module summary:\n";
    for(const auto &i : stats_) {
      os<<"    "<< i.second<< " events with status "<<i.first<<"\n";
    }

    mf::LogInfo("Summary")<<os.str();
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::FilterStatusG4);
