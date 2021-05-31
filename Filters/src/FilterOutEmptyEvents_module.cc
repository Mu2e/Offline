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
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "MCDataProducts/inc/StatusG4.hh"

using namespace std;

namespace mu2e {

  //================================================================
  class FilterOutEmptyEvents : public art::EDFilter {
      int numpassedEvents_ = 0;
  public:
      struct Config {
          fhicl::Atom<string> inputModuleLabel{fhicl::Name("input"), "g4run"};
      };
      
      using Parameters = art::EDFilter::Table<Config>;

  private:
      
      string _inputModuleLabel;
      
  public:
    explicit FilterOutEmptyEvents(Parameters const& config);
    virtual bool filter(art::Event& event) override;
    virtual void endJob() override;
  };
    
  //================================================================
  FilterOutEmptyEvents::FilterOutEmptyEvents(Parameters const& config)
    : art::EDFilter{config},
      _inputModuleLabel(config().inputModuleLabel())
  {}
    
  //================================================================
  bool FilterOutEmptyEvents::filter(art::Event& event) {
      
      bool passed = false;
      
      art::Handle<StatusG4> StatusG4Handle;
      event.getByLabel(_inputModuleLabel,StatusG4Handle);

      //auto ih = event.getHandle(token_);
      
      if (StatusG4Handle.isValid()) {
          passed = true;
          numpassedEvents_ ++;
      }

      return passed;
  }

  //================================================================
  void FilterOutEmptyEvents::endJob() {
    std::ostringstream os;
    os << "FilterOutEmptyEvents_module summary:\n";
    os << "    " << numpassedEvents_ << " events passed \n";

    mf::LogInfo("Summary")<<os.str();
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::FilterOutEmptyEvents);
