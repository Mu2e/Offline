//
//  Study energy desposited in straws 
//
// $Id: StrawEnergy_module.cc,v 1.1 2012/01/20 19:18:19 brownd Exp $
// $Author: 
// $Date: 2012/01/20 19:18:19 $
//
// Original author David Brown
//
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "TTree.h"
#include "TBranch.h"
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

namespace mu2e {

  class StrawEnergy : public art::EDAnalyzer {
  public:

    explicit StrawEnergy(fhicl::ParameterSet const& pset);
    virtual ~StrawEnergy() { }

    virtual void beginJob();
    virtual void endJob();

    // This is called for each event.
    virtual void analyze(const art::Event& e);

  private:
    int _diagLevel;
    std::string _g4ModuleLabel;
    std::string _trackerStepPoints;
    TTree* _estraw;
// branch variables
    Int_t _device, _sector, _layer, _straw;
    Float_t _energy;
    Float_t _time;
    Float_t _x;
    Float_t _rho;
    Float_t _z;
 
  };


  StrawEnergy::StrawEnergy(fhicl::ParameterSet const& pset) :
    // Run time parameters
    _diagLevel(pset.get<int>("diagLevel",0)),
    _g4ModuleLabel(pset.get<string>("g4ModuleLabel","g4run")),
    _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker"))
  {

  }

  void StrawEnergy::beginJob(){
    // Get access to the TFile service.
    art::ServiceHandle<art::TFileService> tfs;
    _estraw=tfs->make<TTree>("estraw","straw energy");
    _estraw->Branch("device",&_device,"device/i");
    _estraw->Branch("sector",&_sector,"sector/i");
    _estraw->Branch("layer",&_layer,"layer/i");
    _estraw->Branch("straw",&_straw,"straw/i");
    _estraw->Branch("energy",&_energy,"energy/f");
    _estraw->Branch("time",&_time,"time/f");
    _estraw->Branch("x",&_x,"x/f");
    _estraw->Branch("rho",&_rho,"rho/f");
    _estraw->Branch("z",&_z,"z/f");
 
  }

  void StrawEnergy::endJob(){
  }

  void StrawEnergy::analyze(const art::Event& event) {
    const Tracker& tracker = getTrackerOrThrow();


    art::Handle<StepPointMCCollection> hits;
    event.getByLabel(_g4ModuleLabel,_trackerStepPoints,hits);
    // Loop over all hits
    for ( size_t i=0; i<hits->size(); ++i ){
      const StepPointMC& hit = (*hits)[i];
      _energy = hit.eDep();
      _time = hit.time();
      // Get the hit information.
      const CLHEP::Hep3Vector& pos = hit.position();
      // Get the straw information:
      const Straw&      straw = tracker.getStraw( hit.strawIndex() );
      _device = straw.id().getDevice();
      _sector = straw.id().getSector();
      _layer = straw.id().getLayer();
      _straw = straw.id().getStraw();
      const CLHEP::Hep3Vector& mid   = straw.getMidPoint();
      const CLHEP::Hep3Vector& w     = straw.getDirection();
// define position using straw
      _x = w.dot(pos-mid);
      _rho = mid.perp();
      _z = mid.z();
      _estraw->Fill();
    }
  }

}

DEFINE_ART_MODULE(mu2e::StrawEnergy);
