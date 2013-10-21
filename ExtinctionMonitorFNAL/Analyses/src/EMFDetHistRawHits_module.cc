// ExtMonFNAL raw hits histograms.
//
// Andrei Gaponenko, 2012

#include <iostream>
#include <string>
#include <vector>

#include "TH2.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "RecoDataProducts/inc/ExtMonFNALRawHitCollection.hh"
#include "ExtinctionMonitorFNAL/Analyses/inc/EMFRawHitHistograms.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "GeometryService/inc/GeomHandle.hh"

namespace mu2e {

  //================================================================
  class EMFDetHistRawHits : public art::EDAnalyzer {
    std::string inputModuleLabel_;
    std::string inputInstanceName_;
    std::string geomModuleLabel_;
    std::string geomInstanceName_;

    // This is a workaround for geometry not being available at beginJob()
    bool booked_;
    EMFRawHitHistograms hh_;

  public:
    explicit EMFDetHistRawHits(const fhicl::ParameterSet& pset);
    virtual void beginRun(const art::Run& run);
    virtual void analyze(const art::Event& event);
  };

  //================================================================
  EMFDetHistRawHits::EMFDetHistRawHits(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset)
    , inputModuleLabel_(pset.get<std::string>("inputModuleLabel"))
    , inputInstanceName_(pset.get<std::string>("inputInstanceName", ""))
    , geomModuleLabel_(pset.get<std::string>("geomModuleLabel"))
    , geomInstanceName_(pset.get<std::string>("geomInstanceName", ""))
    , booked_(false)
    , hh_(pset)
  {}

  //================================================================
  void EMFDetHistRawHits::beginRun(const art::Run& run) {
    // This is a workaround for geometry not being available at beginJob()
    if(!booked_) {
      booked_ = true;
      if(!geomModuleLabel_.empty()) {
        art::Handle<ExtMonFNAL::ExtMon> extmon;
        run.getByLabel(geomModuleLabel_, geomInstanceName_, extmon);
        hh_.book(*extmon);
      }
      else {
        GeomHandle<ExtMonFNAL::ExtMon> extmon;
        hh_.book(*extmon);
      }
    }
  }

  //================================================================
  void EMFDetHistRawHits::analyze(const art::Event& event) {
    art::Handle<ExtMonFNALRawHitCollection> ih;
    event.getByLabel(inputModuleLabel_, inputInstanceName_, ih);
    hh_.fill(*ih);
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::EMFDetHistRawHits);
