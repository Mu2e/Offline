//
// example Plugin to read Tracker PlaneSupport Detectors data and create ntuples
//
//  $Id: ReadTrackerDSDetectors_module.cc,v 1.2 2013/10/21 20:44:04 genser Exp $
//  $Author: genser $
//  $Date: 2013/10/21 20:44:04 $
//
// Original author KLG
//

//#include <cmath>
#include <iostream>
#include <string>
#include <iomanip>


// CLHEP includes
#include "CLHEP/Units/SystemOfUnits.h"

// Mu2e includes
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

// root includes
#include "TH1F.h"
#include "TNtuple.h"
#include "TTree.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

using namespace std;

using CLHEP::Hep3Vector;
using CLHEP::keV;

namespace mu2e {

  class ReadTrackerDSDetectors : public art::EDAnalyzer {
  public:

    explicit ReadTrackerDSDetectors(fhicl::ParameterSet const& pset);
    virtual ~ReadTrackerDSDetectors() {}

    virtual void beginJob();
    void analyze(const art::Event& e);

  private:

    // debug output llevel
    int _diagLevel;

   // Module label of the module that made the hits.
    std::string _hitMakerModuleLabel;

    // Control printed output.
    int _nAnalyzed;
    int _maxFullPrint;

    // Pointers to histograms, ntuples.

    TH1F*    _hNsdDet;
    TH1F*    _hNsdDetH;

    TNtuple* _nttsdd;

    // Name of the SDD StepPoint collection
    std::string  _sddStepPoints;

  };

  ReadTrackerDSDetectors::ReadTrackerDSDetectors(fhicl::ParameterSet const& pset) : 
    art::EDAnalyzer(pset),
    // Run time parameters
    _diagLevel(pset.get<int>("diagLevel",0)),
    _hitMakerModuleLabel(pset.get<string>("hitMakerModuleLabel","g4run")),
    _nAnalyzed(0),
    _maxFullPrint(pset.get<int>("maxFullPrint",5)),
    _hNsdDet(0),
    _hNsdDetH(0),
    _nttsdd(0),
    _sddStepPoints(pset.get<string>("sddStepPoints","trackerDS"))
  {}

 

  void ReadTrackerDSDetectors::beginJob(){

    // Get access to the TFile service.

    art::ServiceHandle<art::TFileService> tfs;

    _hNsdDet  = tfs->make<TH1F>( "hNsdDet ", 
                                 "Number/ID of the detector which was hit",
                                 50,  0., 50. );

    _hNsdDetH = tfs->make<TH1F>( "hNsdDetH", 
                                 "Number of hits in PlaneSupport Detectors",
                                 50,  0., 50. );

    // Create an ntuple.
    _nttsdd = tfs->make<TNtuple>( "nttsdd", 
                                  "Plane Support Detectors ntuple",
                                  "evt:trk:sid:pdg:time:x:y:z:px:py:pz:"
                                  "g4bl_time");
  }

  void ReadTrackerDSDetectors::analyze(const art::Event& event) {

    ++_nAnalyzed;

    if (_diagLevel>1 ) {
      cout << "ReadTrackerDSDetectors::" << __func__ 
           << setw(4) << " called for event "  
           << event.id().event()
           << " hitMakerModuleLabel "
           << _hitMakerModuleLabel 
           << " sddStepPoints " 
           << _sddStepPoints
           << endl;
    }

    art::ServiceHandle<GeometryService> geom;

    if (!geom->hasElement<Tracker>()) 
      {
        mf::LogError("Geom")
          << "Skipping ReadTrackerDSDetectors::analyze due to lack of tracker\n";
        return;
      }

    // Access detector geometry information

    // Get a handle to the hits created by G4.
    art::Handle<StepPointMCCollection> hitsHandle;

    event.getByLabel(_hitMakerModuleLabel, _sddStepPoints, hitsHandle);

    if (!hitsHandle.isValid() ) {

      mf::LogError("Hits")
        << " Skipping ReadTrackerDSDetectors::analyze due to problems with hits\n";
      return;
    }

    StepPointMCCollection const& hits = *hitsHandle;

    unsigned int const nhits = hits.size();

    if (_diagLevel>0 && nhits>0 ) {
      cout << "ReadTrackerDSDetectors::" << __func__ 
           << " Number of SDD hits: " <<setw(4) 
           << nhits
           << endl;
    }

    art::Handle<SimParticleCollection> simParticles;
    event.getByLabel(_hitMakerModuleLabel, simParticles);
    bool haveSimPart = simParticles.isValid();
    if ( haveSimPart ) haveSimPart = !(simParticles->empty());

    // Fill histograms & ntuple

    // Fill histogram with number of hits per event.
    _hNsdDetH->Fill(nhits);

    float nt[_nttsdd->GetNvar()];

    // Loop over all SDD hits.
    for ( size_t i=0; i<nhits; ++i ){

      // Alias, used for readability.
      const StepPointMC& hit = (hits)[i];

      // Get the hit information.

      int did = hit.volumeId();

      // Fill histogram with the detector numbers
      _hNsdDet->Fill(did);

      const CLHEP::Hep3Vector& pos = hit.position();
      const CLHEP::Hep3Vector& mom = hit.momentum();

      // Get track info
      SimParticleCollection::key_type trackId = hit.trackId();
      int pdgId = 0;
      if ( haveSimPart ){
        if( !simParticles->has(trackId) ) {
          pdgId = 0;
        } else {
          SimParticle const& sim = simParticles->at(trackId);
          pdgId = sim.pdgId();
        }
      }

      // Fill the ntuple.
      nt[0]  = event.id().event();
      nt[1]  = trackId.asInt();
      nt[2]  = did;
      nt[3]  = pdgId;
      nt[4]  = hit.time();
      nt[5]  = pos.x();
      nt[6]  = pos.y();
      nt[7]  = pos.z();
      nt[8]  = mom.x();
      nt[9]  = mom.y();
      nt[10] = mom.z();
      nt[11] = hit.properTime();

      _nttsdd->Fill(nt);

      if ( _nAnalyzed < _maxFullPrint){
        cout <<  "ReadTrackerDSDetectors::" << __func__ 
             << ": SDD hit: "
             << setw(8) << event.id().event() << " | "
             << setw(4) << hit.volumeId()     << " | "
             << setw(6) << pdgId              << " | "
             << setw(8) << hit.time()         << " | "  
             << setw(8) << mom.mag()          << " | "
             << pos               
             << endl;

      }

    } // end loop over hits.

  }

}  // end namespace mu2e

using mu2e::ReadTrackerDSDetectors;
DEFINE_ART_MODULE(ReadTrackerDSDetectors);
