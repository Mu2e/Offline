//
// Plugin to test that I can read back the persistent data
// from crude straw hits.  Also tests:
//   - CrudeStrawHitCollection
//   - the mechanisms to look back at the precursor StepPointMC objects.
//
// $Id: MCSH_Test_plugin.cc,v 1.9 2010/09/30 02:10:26 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/09/30 02:10:26 $
//
// Original author Rob Kutschke
//

// C++ includes.
#include <iostream>
#include <iomanip>
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
#include "LTrackerGeom/inc/LTracker.hh"
#include "LTrackerGeom/inc/CrudeStrawHitCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "Mu2eUtilities/inc/resolveTransients.hh"


using namespace std;

namespace mu2e {

  //--------------------------------------------------------------------
  //
  // 
  class MCSH_Test : public edm::EDAnalyzer {
  public:
    explicit MCSH_Test(edm::ParameterSet const& pset):
      _diagLevel(pset.getUntrackedParameter<int>("diagLevel",0)),
      _maxFullPrint(pset.getUntrackedParameter<int>("maxFullPrint",5)),
      _hDriftDist(0){
    }
    virtual ~MCSH_Test() { }

    virtual void beginJob(edm::EventSetup const&);

    void analyze( edm::Event const& e, edm::EventSetup const&);

  private:
    
    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Some diagnostic histograms.
    TH1F* _hDriftDist;

  };

  void MCSH_Test::beginJob(edm::EventSetup const& ){

    cout << "Diaglevel: " 
         << _diagLevel << " "
         << _maxFullPrint 
         << endl;

    // Create histograms if diagnostics are enabled.
    if ( _diagLevel > 0 ){

      edm::Service<edm::TFileService> tfs;

      _hDriftDist = tfs->make<TH1F>( "hDriftDist", "True Drift Distance;(mm)", 100,  0.,   3.  );

    }

  }

  void
  MCSH_Test::analyze(edm::Event const& evt, edm::EventSetup const&) {

    static int ncalls(0);
    ++ncalls;

    // Instance name of the module that created the hits of interest;
    static const string creatorName("makeCSH");

    // Geometry info for the LTracker.
    GeomHandle<LTracker> ltracker;

    // Get the persistent data about the CrudeStrawHits.
    edm::Handle<CrudeStrawHitPData> pdataHandle;
    evt.getByLabel(creatorName,pdataHandle);
    CrudeStrawHitPData const* pdata = pdataHandle.product();

    if ( ncalls < _maxFullPrint && _diagLevel > 2 ){
      for ( std::size_t i=0; i<pdata->size(); ++i){
        cout << "Before pdata: " 
             << evt.id().event() <<  " " 
             << i << " " 
             << pdata->at(i).stepPointMCsValid() << " "
             << pdata->at(i).getStepPointMCs(true).size() << " "
             << pdata->at(i).getStepPointMCs(true).capacity() << " "
             << endl;
      }
    }

    // Form a fully functional collection of these hits.
    //    CrudeStrawHitCollection crudeHits(evt, pdata);
    //crudeHits.resolveTransients(evt);

    resolveTransients<CrudeStrawHitPData>( *pdata, evt); 

    if ( ncalls < _maxFullPrint && _diagLevel > 2 ){
      for ( std::size_t i=0; i<pdata->size(); ++i){
        cout << "After pdata: " 
             << evt.id().event() <<  " " 
             << i << " " 
             << pdata->at(i).stepPointMCsValid() << " "
             << pdata->at(i).getStepPointMCs().size() << " "
             << pdata->at(i).getStepPointMCs().capacity() << " "
             << endl;
      }
    }


    for ( vector<CrudeStrawHit>::size_type i=0; 
          i<pdata->size(); 
          ++i){

      // Aliases for readability.
      CrudeStrawHit const&      hit(pdata->at(i));
      Straw const&              straw(ltracker->getStraw(hit.strawIndex));

      // The list of nearest neighbours of this straw.
      //vector<StrawIndex> const& nearest(straw.nearestNeighboursByIndex());

      // Fill diagnostic histogram.
      if ( _diagLevel > 0){
        _hDriftDist->Fill(hit.trueDriftDistance);
      }
      
      // Diagnostic printout.
      if ( ncalls < _maxFullPrint && _diagLevel > 2 ){
        
        // Print list of nearest neighbours that are hit.
        cout << "Hit neighbours : " 
             << setw(4) << i <<  " : "
             << setw(4) << straw.Id() << " : ";

        /*
        for ( vector<StrawIndex>::size_type j=0;
              j<nearest.size(); 
              ++j ){

          if ( crudeHits.hasHitByStrawIndex(nearest[j]) ){
            Straw const& s = ltracker->getStraw(nearest[j]);
            cout << " " << setw(4) << s.Id();
          }

        }
        */
        cout << endl;

        // Get pointers back to precursors of this hit.
        vector<StepPointMC const*> v;
        pdata->at(i).getStepPointMC(evt, v);
        
        cout << "Roundtrip Id Check for straw: "
             << hit.strawIndex << " | #points: "
             << v.size() << " | points: ";
        for ( vector<StepPointMC const*>::size_type i=0;
              i<v.size(); ++i ){
          cout << " " << v[i]->strawIndex();
        }
        cout << endl;

        // Exercise two other versions of the DPIndex resolver.
        // Should move this testing to ToyDP/test.
        if ( hit.precursorIndices.size() > 0 ) {
          StepPointMC const * p = resolveDPIndex<StepPointMCCollection>( evt, hit.precursorIndices[0]);

          vector<StepPointMC const*> v2;
          edm::ProductID const& id = hit.precursorIndices[0].id;
          vector<int> offsets(1,hit.precursorIndices[0].index);
          resolveDPIndices<StepPointMCCollection>( evt, id, offsets, v2);

          cout << "First precursor check: "
               << hit.strawIndex  << " : "
               << p->strawIndex() << "   "
               << v2[0]->strawIndex()
               << endl;
        }

      } // end of ncalls < _maxFullPrint

    } // end of loop over hits.
    
  } // end of ::analyze.
  
}


using mu2e::MCSH_Test;
DEFINE_FWK_MODULE(MCSH_Test);
