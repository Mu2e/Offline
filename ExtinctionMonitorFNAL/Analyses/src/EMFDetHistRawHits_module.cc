// ExtMonFNAL raw hits histograms.
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

#include "RecoDataProducts/inc/ExtMonFNALRawHitCollection.hh"
#include "ExtinctionMonitorFNAL/Analyses/inc/EMFRawHitHistograms.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "GeometryService/inc/GeomHandle.hh"

namespace mu2e {

  //================================================================
  class EMFDetHistRawHits : public art::EDAnalyzer {
    std::string _inModuleLabel;
    std::string _inInstanceName;

    EMFRawHitHistograms hh_;

  public:
    explicit EMFDetHistRawHits(const fhicl::ParameterSet& pset);
    virtual void beginRun(const art::Run& run);
    virtual void analyze(const art::Event& event);
  };

  //================================================================
  EMFDetHistRawHits::EMFDetHistRawHits(const fhicl::ParameterSet& pset)
    : _inModuleLabel(pset.get<std::string>("inputModuleLabel"))
    , _inInstanceName(pset.get<std::string>("inputInstanceName"))
  {}

  //================================================================
  void EMFDetHistRawHits::beginRun(const art::Run&) {
    GeomHandle<ExtMonFNAL::ExtMon> extmon;
    hh_.book(*extmon);
  }

  //================================================================
  void EMFDetHistRawHits::analyze(const art::Event& event) {
    art::Handle<ExtMonFNALRawHitCollection> ih;
    event.getByLabel(_inModuleLabel, _inInstanceName, ih);
    hh_.fill(*ih);
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::EMFDetHistRawHits);
