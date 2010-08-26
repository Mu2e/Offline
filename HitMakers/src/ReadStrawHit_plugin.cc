//
// Plugin to test that I can read back the persistent data about straw hits.  
// Also tests the mechanisms to look back at the precursor StepPointMC objects.
//
// $Id: ReadStrawHit_plugin.cc,v 1.3 2010/08/26 19:16:45 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/08/26 19:16:45 $
//
// Original author Rob Kutschke. Updated by Ivan Logashenko.
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

// Framework includes.
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Root includes.
#include "TFile.h"
#include "TH1F.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "ToyDP/inc/StrawHitCollection.hh"
#include "ToyDP/inc/StrawHitMCTruthCollection.hh"
#include "ToyDP/inc/StrawHitMCPtrCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "Mu2eUtilities/inc/resolveTransients.hh"

// Other includes.
#include "CLHEP/Random/RandGauss.h"


using namespace std;
using edm::Event;

namespace mu2e {

  //--------------------------------------------------------------------
  //
  // 
  class ReadStrawHit : public edm::EDAnalyzer {
  public:
    explicit ReadStrawHit(edm::ParameterSet const& pset):
      _diagLevel(pset.getUntrackedParameter<int>("diagLevel",0)),
      _maxFullPrint(pset.getUntrackedParameter<int>("maxFullPrint",5)),
      _makerModuleLabel(pset.getParameter<std::string>("makerModuleLabel")),
      _hHitTime(0),
      _hHitDeltaTime(0),
      _hHitEnergy(0),
      _hNHits(0),
      _hNHitsPerWire(0),
      _hDriftTime(0),
      _hDriftDistance(0),
      _hDistanceToMid(0),
      _hNG4Steps(0),
      _hT0(0),
      _hG4StepLength(0),
      _hG4StepEdep(0)
    {
    }
    virtual ~ReadStrawHit() { }

    virtual void beginJob(edm::EventSetup const&);

    void analyze( edm::Event const& e, edm::EventSetup const&);

  private:
    
    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Label of the module that made the hits.
    std::string _makerModuleLabel;

    // Some diagnostic histograms.
    TH1F* _hHitTime;
    TH1F* _hHitDeltaTime;
    TH1F* _hHitEnergy;
    TH1F* _hNHits;
    TH1F* _hNHitsPerWire;
    TH1F* _hDriftTime;
    TH1F* _hDriftDistance;
    TH1F* _hDistanceToMid;
    TH1F* _hNG4Steps;
    TH1F* _hT0;
    TH1F* _hG4StepLength;
    TH1F* _hG4StepEdep;

  };

  void ReadStrawHit::beginJob(edm::EventSetup const& ){

    cout << "Diaglevel: " 
         << _diagLevel << " "
         << _maxFullPrint 
         << endl;

    edm::Service<edm::TFileService> tfs;

    _hHitTime      = tfs->make<TH1F>( "hHitTime",      "Hit Time (ns)", 200, 0., 2000. );
    _hHitDeltaTime = tfs->make<TH1F>( "hHitDeltaTime", "Hit Delta Time (ns)", 80, -20.0, 20. );
    _hHitEnergy    = tfs->make<TH1F>( "hHitEnergy",    "Hit Energy (keV)", 100, 0., 100. );
    _hNHits        = tfs->make<TH1F>( "hNHits",        "Number of straw hits", 500, 0., 500. );
    _hNHitsPerWire = tfs->make<TH1F>( "hNHitsPerWire", "Number of hits per straw", 10, 0., 10. );
    _hDriftTime    = tfs->make<TH1F>( "hDriftTime",    "Drift time, ns", 100, 0., 100. );
    _hDriftDistance= tfs->make<TH1F>( "hDriftDistance","Drift Distance, mm", 100, 0., 3. );
    _hDistanceToMid= tfs->make<TH1F>( "hDistanceToMid","Distance to wire center, mm", 160, -1600., 1600. );
    _hNG4Steps     = tfs->make<TH1F>( "hNG4Steps",     "Number of G4Steps per hit", 100, 0., 100. );
    _hT0           = tfs->make<TH1F>( "hT0",           "T0, ns", 100, -50., 50. );
    _hG4StepLength = tfs->make<TH1F>( "hG4StepLength", "Length of G4Steps, mm", 100, 0., 10. );
    _hG4StepEdep   = tfs->make<TH1F>( "hG4StepEdep",   "Energy deposition of G4Steps, keV", 100, 0., 10. );

  }

  void
  ReadStrawHit::analyze(edm::Event const& evt, edm::EventSetup const&) {

    static int ncalls(0);
    ++ncalls;

    /*
    // Print the content of current event
    std::vector<edm::Provenance const*> edata;
    evt.getAllProvenance(edata);
    cout << "Event info: " 
    << evt.id().event() <<  " " 
    << " found " << edata.size() << " objects " 
    << endl;
    for( int i=0; i<edata.size(); i++ ) cout << *(edata[i]) << endl;
    */

    // Geometry info for the LTracker.
    // Get a reference to one of the L or T trackers.
    // Throw exception if not successful.

    const Tracker& tracker = getTrackerOrThrow();

    // Get the persistent data about the StrawHits.

    edm::Handle<StrawHitCollection> pdataHandle;
    evt.getByLabel(_makerModuleLabel,pdataHandle);
    StrawHitCollection const* hits = pdataHandle.product();

    // Get the persistent data about the StrawHitsMCTruth.

    edm::Handle<StrawHitMCTruthCollection> truthHandle;
    evt.getByLabel(_makerModuleLabel,truthHandle);
    StrawHitMCTruthCollection const* hits_truth = truthHandle.product();

    // Get the persistent data about pointers to StepPointMCs

    edm::Handle<StrawHitMCPtrCollection> mcptrHandle;
    evt.getByLabel(_makerModuleLabel,mcptrHandle);
    StrawHitMCPtrCollection const* hits_mcptr = mcptrHandle.product();

    // Get the persistent data about the StepPointMCs. More correct implementation
    // should look for product ids in StrawHitMCPtrCollection, rather than 
    // use producer name directly ("g4run"). 

    edm::Handle<StepPointMCCollection> mchitsHandle;
    evt.getByLabel("g4run",mchitsHandle);
    StepPointMCCollection const* mchits = mchitsHandle.product();

    // Fill histograms

    _hNHits->Fill(hits->size());

    std::map<StrawIndex,int> nhperwire;

    for ( int i=0; i<hits->size(); ++i ) {

      // Access data
      StrawHit        const&      hit(hits->at(i));
      StrawHitMCTruth const&    truth(hits_truth->at(i));
      StrawHitMCPtr   const&    mcptr(hits_mcptr->at(i));
      
      // Fill per-event histograms
      if( i==0 ) {
        _hT0->Fill(truth.t0());
      }

      // Use data from hits
      _hHitTime->Fill(hit.time());
      _hHitDeltaTime->Fill(hit.dt());
      _hHitEnergy->Fill(hit.energyDep()*1000.0);

      // Use MC truth data
      _hDriftTime->Fill(truth.driftTime());
      _hDriftDistance->Fill(truth.driftDistance());
      _hDistanceToMid->Fill(truth.distanceToMid());

      // Use data from G4 hits
      _hNG4Steps->Fill(mcptr.size());
      for( int j=0; j<mcptr.size(); ++j ) {
        StepPointMC const& mchit = (*mchits)[mcptr[j].index];
        _hG4StepLength->Fill(mchit.stepLength());
        _hG4StepEdep->Fill(mchit.eDep()*1000.0);
      }

      // Calculate number of hits per wire
      nhperwire[hit.strawIndex()]++;

    }

    for( std::map<StrawIndex,int>::iterator it=nhperwire.begin(); it!= nhperwire.end(); ++it ) {
      _hNHitsPerWire->Fill(it->second);
    }

  } // end of ::analyze.
  
}


using mu2e::ReadStrawHit;
DEFINE_FWK_MODULE(ReadStrawHit);
