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
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"

#include "RecoDataProducts/inc/ExtMonFNALRecoClusterCollection.hh"
#include "ExtinctionMonitorFNAL/Analyses/inc/EMFRecoClusterHistograms.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "GeometryService/inc/GeomHandle.hh"

namespace mu2e {

  //================================================================
  class EMFDetHistRecoClusters : public art::EDAnalyzer {
    std::string inModuleLabel_;
    std::string inInstanceName_;
    std::string geomModuleLabel_;
    std::string geomInstanceName_;

    // This is a workaround for geometry not being available at beginJob()
    bool booked_;
    EMFRecoClusterHistograms ch_;

  public:
    explicit EMFDetHistRecoClusters(const fhicl::ParameterSet& pset);
    virtual void beginRun(const art::Run& run);
    virtual void analyze(const art::Event& event);
  };

  //================================================================
  EMFDetHistRecoClusters::EMFDetHistRecoClusters(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset)
    , inModuleLabel_(pset.get<std::string>("inputModuleLabel"))
    , inInstanceName_(pset.get<std::string>("inputInstanceName", ""))
    , geomModuleLabel_(pset.get<std::string>("geomModuleLabel"))
    , geomInstanceName_(pset.get<std::string>("geomInstanceName", ""))
    , booked_(false)
  {}

  //================================================================
  void EMFDetHistRecoClusters::beginRun(const art::Run& run) {
    // This is a workaround for geometry not being available at beginJob()
    if(!booked_) {
      booked_ = true;
      if(!geomModuleLabel_.empty()) {
        art::Handle<ExtMonFNAL::ExtMon> emf;
        run.getByLabel(geomModuleLabel_, geomInstanceName_, emf);
        ch_.book(*emf);
      }
      else {
        GeomHandle<ExtMonFNAL::ExtMon> emf;
        ch_.book(*emf);
      }
    }
  }

  //================================================================
  void EMFDetHistRecoClusters::analyze(const art::Event& event) {
    art::Handle<ExtMonFNALRecoClusterCollection> ih;
    event.getByLabel(inModuleLabel_, inInstanceName_, ih);
    ch_.fill(*ih);
  }

  //================================================================
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::EMFDetHistRecoClusters);
