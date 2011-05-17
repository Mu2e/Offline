//
// Plugin to test that I can read back the persistent data
// from crude straw hits.  Also tests:
//   - CrudeStrawHitCollection
//   - the mechanisms to look back at the precursor StepPointMC objects.
//
// $Id: MCSH_Test_module.cc,v 1.1 2011/05/17 16:30:14 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 16:30:14 $
//
// Original author Rob Kutschke
//

// C++ includes.
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

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
  class MCSH_Test : public art::EDAnalyzer {
  public:
    explicit MCSH_Test(fhicl::ParameterSet const& pset):
      _diagLevel(pset.get<int>("diagLevel",0)),
      _maxFullPrint(pset.get<int>("maxFullPrint",5)),
      _hDriftDist(0){
    }
    virtual ~MCSH_Test() { }

    virtual void beginJob(art::EventSetup const&);

    void analyze( art::Event const& e, art::EventSetup const&);

  private:
    
    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Some diagnostic histograms.
    TH1F* _hDriftDist;

  };

  void MCSH_Test::beginJob(art::EventSetup const& ){

    cout << "Diaglevel: " 
         << _diagLevel << " "
         << _maxFullPrint 
         << endl;

    // Create histograms if diagnostics are enabled.
    if ( _diagLevel > 0 ){

      art::ServiceHandle<art::TFileService> tfs;

      _hDriftDist = tfs->make<TH1F>( "hDriftDist", "True Drift Distance;(mm)", 100,  0.,   3.  );

    }

  }

  void
  MCSH_Test::analyze(art::Event const& evt, art::EventSetup const&) {

    static int ncalls(0);
    ++ncalls;

    // Instance name of the module that created the hits of interest;
    static const string creatorName("makeCSH");

    // Geometry info for the LTracker.
    GeomHandle<LTracker> ltracker;

    // Get the persistent data about the CrudeStrawHits.
    art::Handle<CrudeStrawHitPData> pdataHandle;
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
          art::ProductID const& id = hit.precursorIndices[0].id;
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
DEFINE_ART_MODULE(MCSH_Test);
