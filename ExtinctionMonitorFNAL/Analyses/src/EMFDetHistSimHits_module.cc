// ExtMonFNAL cluster histograms.
//
// Andrei Gaponenko, 2012

#include <iostream>
#include <string>
#include <vector>

#include "TH2.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "MCDataProducts/inc/ExtMonFNALSimHitCollection.hh"
#include "ExtinctionMonitorFNAL/Analyses/inc/EMFSimHitHistograms.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "GeometryService/inc/GeomHandle.hh"

namespace mu2e {

  //================================================================
  class EMFDetHistSimHits : public art::EDAnalyzer {
    std::string _inModuleLabel;
    std::string _inInstanceName;

    EMFSimHitHistograms ch_;

  public:
    explicit EMFDetHistSimHits(const fhicl::ParameterSet& pset);
    virtual void beginRun(const art::Run& run);
    virtual void analyze(const art::Event& event);
  };

  //================================================================
  EMFDetHistSimHits::EMFDetHistSimHits(const fhicl::ParameterSet& pset)
    : _inModuleLabel(pset.get<std::string>("inputModuleLabel"))
    , _inInstanceName(pset.get<std::string>("inputInstanceName"))
  {}

  //================================================================
  void EMFDetHistSimHits::beginRun(const art::Run&) {
    GeomHandle<ExtMonFNAL::ExtMon> extmon;
    ch_.book(*extmon);
  }

  //================================================================
  void EMFDetHistSimHits::analyze(const art::Event& event) {
    art::Handle<ExtMonFNALSimHitCollection> ih;
    event.getByLabel(_inModuleLabel, _inInstanceName, ih);
    ch_.fill(*ih);
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::EMFDetHistSimHits);
