// Pass events with StatusG4 <= maxAcceptedStatus.
// See MCDataProducts/inc/StatusG4.hh for the meaning of different status values.
// By default maxAcceptedStatus=0 and only perfect events are passed.
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
  class FilterCRYOut : public art::EDFilter {
      art::InputTag input_;
      int numpassedEvents_ = 0;
  public:
    explicit FilterCRYOut(const fhicl::ParameterSet& pset);
    virtual bool filter(art::Event& event) override;
    virtual void endJob() override;
  };

  //================================================================
  FilterCRYOut::FilterCRYOut(const fhicl::ParameterSet& pset)
    : input_(pset.get<std::string>("input","g4un"))
  {}

  //================================================================
  bool FilterCRYOut::filter(art::Event& event) {
      
      bool passed = false;
      
      art::Handle<StatusG4> StatusG4Handle;
      event.getByLabel(input_,StatusG4Handle);

      //auto ih = event.getHandle(token_);
      
      if (StatusG4Handle.isValid()) {
          passed = true;
          numpassedEvents_ ++;
      }

      return passed;
  }

  //================================================================
  void FilterCRYOut::endJob() {
    std::ostringstream os;
    os << "FilterCRYOut_module summary:\n";
    os << "    " << numpassedEvents_ << " events passed \n";

    mf::LogInfo("Summary")<<os.str();
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::FilterCRYOut);
