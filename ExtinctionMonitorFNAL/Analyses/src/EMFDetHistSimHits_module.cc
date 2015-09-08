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

#include "MCDataProducts/inc/ExtMonFNALSimHitCollection.hh"
#include "ExtinctionMonitorFNAL/Analyses/inc/EMFSimHitHistograms.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "GeometryService/inc/GeomHandle.hh"

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    class EMFDetHistSimHits : public art::EDAnalyzer {
      std::string _inModuleLabel;
      std::string _inInstanceName;
      std::string geomModuleLabel_;
      std::string geomInstanceName_;

      // This is a workaround for geometry not being available at beginJob()
      bool booked_;
      EMFSimHitHistograms ch_;


    public:
      explicit EMFDetHistSimHits(const fhicl::ParameterSet& pset);
      virtual void beginRun(const art::Run& run);
      virtual void analyze(const art::Event& event);
    };

    //================================================================
    EMFDetHistSimHits::EMFDetHistSimHits(const fhicl::ParameterSet& pset)
      : art::EDAnalyzer(pset)
      , _inModuleLabel(pset.get<std::string>("inputModuleLabel"))
      , _inInstanceName(pset.get<std::string>("inputInstanceName", ""))
      , geomModuleLabel_(pset.get<std::string>("geomModuleLabel"))
      , geomInstanceName_(pset.get<std::string>("geomInstanceName", ""))
      , booked_(false)
    {}

    //================================================================
    void EMFDetHistSimHits::beginRun(const art::Run& run) {
      // This is a workaround for geometry not being available at beginJob()
      if(!booked_) {
        booked_ = true;
        if(!geomModuleLabel_.empty()) {
          art::Handle<ExtMonFNAL::ExtMon> extmon;
          run.getByLabel(geomModuleLabel_, geomInstanceName_, extmon);
          ch_.book(*extmon);
        }
        else {
          GeomHandle<ExtMonFNAL::ExtMon> extmon;
          ch_.book(*extmon);
        }
      }
    }

    //================================================================
    void EMFDetHistSimHits::analyze(const art::Event& event) {
      GeomHandle<ExtMonFNAL::ExtMon> extmon;
      art::Handle<ExtMonFNALSimHitCollection> ih;
      event.getByLabel(_inModuleLabel, _inInstanceName, ih);
      ch_.fill(*extmon, *ih);
    }

    //================================================================
  } // namespace ExtMonFNAL
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::EMFDetHistSimHits);
