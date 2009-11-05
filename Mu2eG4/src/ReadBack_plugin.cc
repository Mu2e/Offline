//
// An EDAnalyzer Module that reads back the hits created by G4 and makes histograms.
//
// $Id: ReadBack_plugin.cc,v 1.5 2009/11/05 14:38:40 rhbob Exp $
// $Author: rhbob $ 
// $Date: 2009/11/05 14:38:40 $
//
// Original author Rob Kutschke
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
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

// Root includes.
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"

// Mu2e includes.
#include "LTrackerGeom/inc/LTracker.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {

  //--------------------------------------------------------------------
  //
  // 
  class ReadBack : public edm::EDAnalyzer {
  public:
    explicit ReadBack(edm::ParameterSet const& pset) : 
      _maxFullPrint(pset.getUntrackedParameter<int>("maxFullPrint",5)),
      _nAnalyzed(0),
      _hRadius(0),
      _hTime(0),
      _hMultiplicity(0),
      _hDriftDist(0),
      _messageCategory("ToyHitInfo"){
    }
    virtual ~ReadBack() { }

    virtual void beginJob(edm::EventSetup const&);
    virtual void endJob();

    virtual void beginRun(edm::Run const &r, 
			  edm::EventSetup const& eSetup );

    virtual void beginLuminosityBlock(edm::LuminosityBlock const& lblock, 
				      edm::EventSetup const&);
 
    // This is called for each event.
    void analyze(const edm::Event& e, edm::EventSetup const&);

  private:

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Number of events analyzed.
    int _nAnalyzed;

    // Pointers to histograms to be filled.
    TH1F* _hRadius;
    TH1F* _hTime;
    TH1F* _hMultiplicity;
    TH1F* _hDriftDist;
    TH1F* _hxHit;
    TH1F* _hyHit;
    TH1F* _hzHit;
    TH1F* _hHitNeighbours;
    TH1F* _hCheckPointRadius;

    TNtuple* _ntup;

    // A category for the error logger.
    const std::string _messageCategory;

    // A helper function.
    int countHitNeighbours( Straw const& straw, 
			    edm::Handle<StepPointMCCollection>& hits );


  };

  void ReadBack::beginJob(edm::EventSetup const& ){

    // Get access to the TFile service.
    edm::Service<edm::TFileService> tfs;

    // Create some 1D histograms.
    _hRadius       = tfs->make<TH1F>( "hRadius", "Radius of Hits;(mm)",          100,  0., 1000. );
    _hTime         = tfs->make<TH1F>( "hTime", "Pulse Height;(ns)",              100,  0., 2000. );
    _hMultiplicity = tfs->make<TH1F>( "hMultiplicity", "Hits per Event",         100,  0.,  100. );
    _hDriftDist    = tfs->make<TH1F>( "hDriftDist", "Crude Drift Distance;(mm)", 100,  0.,   3.  );
    _hxHit         = tfs->make<TH1F>( "hxHit",  "X of Hit;(mm)",                 
				      100,  -1000.,  1000. );
    _hyHit         = tfs->make<TH1F>( "hyHit",  "Y of Hit;(mm)",                 
				      100,  -1000.,  1000. );
    _hzHit         = tfs->make<TH1F>( "hzHit",  "Z of Hit;(mm)",                 
				      100,  -1400.,  1400. );

    _hHitNeighbours    = tfs->make<TH1F>( "hHitNeighbours",  "Number of hit neighbours",
					  10, 0., 10. );

    _hCheckPointRadius = tfs->make<TH1F>( "hCheckPointRadius",  "Radius of Reference point; (mm)",
					  100, 2.25, 2.75 );

    // Create an ntuple.
    _ntup           = tfs->make<TNtuple>( "ntup", "Hit ntuple", 
					  "evt:trk:sid:hx:hy:hz:wx:wy:wz:dca:time:dev:sec");

  }

  void ReadBack::endJob(){
  }


  void ReadBack::beginRun(edm::Run const& run,
                                 edm::EventSetup const& eSetup ){
  }

  void ReadBack::beginLuminosityBlock(edm::LuminosityBlock const& lblock,
					     edm::EventSetup const&){
  }


  void
  ReadBack::analyze(const edm::Event& evt, edm::EventSetup const&) {

    static int ncalls(0);
    ++ncalls;

    // Master geometry for the LTracker.
    GeomHandle<LTracker> ltracker;

    // Maintain a counter for number of events seen.
    ++_nAnalyzed;

    // Instance name of the module that created the hits of interest;
    static const string creatorName("g4run");

    // Ask the event to give us a "handle" to the requested hits.
    edm::Handle<StepPointMCCollection> hits;
    evt.getByLabel(creatorName,hits);
    
    // Fill histogram with number of hits per event.
    _hMultiplicity->Fill(hits->size());

    // ntuple buffer.
    float nt[13];

    // Loop over all hits.
    int n(0);
    
    for ( StepPointMCCollection::const_iterator 
	    i = hits->begin(), e = hits->end();
	  i!=e; ++i){
      
      // Aliases, used for readability.
      const StepPointMC& hit = *i;
      const Hep3Vector& pos = hit.position();
      const Hep3Vector& mom = hit.momentum();
      
      // Get the straw information.
      Straw const& straw = ltracker->getStraw( hit.strawIndex() );
      Hep3Vector mid = straw.getMidPoint();
      Hep3Vector w   = straw.getDirection();

      // Count how many nearest neighbours are also hit.
      int nNeighbours = countHitNeighbours( straw, hits );

      // Compute an estimate of the drift distance.
      TwoLinePCA pca( mid, w, pos, mom);

      // Check that the radius of the reference point in the local
      // coordinates of the straw.  Should be 2.5 mm.
      double s = w.dot(pos-mid);
      Hep3Vector point = pos - (mid + s*w);

      // I don't understand the distribution of the time variable.
      // I want it to be the time from the start of the spill.
      // It appears to be the time since start of tracking.

      // Fill some histograms
      _hRadius->Fill(pos.perp());
      _hTime->Fill(hit.time());
      _hHitNeighbours->Fill(nNeighbours);
      _hCheckPointRadius->Fill(point.mag());

      _hxHit->Fill(pos.x());
      _hyHit->Fill(pos.y());
      _hzHit->Fill(pos.z());

      _hDriftDist->Fill(pca.dca());

      // Fill the ntuple.
      nt[0]  = evt.id().event();
      nt[1]  = hit.trackId();
      nt[2]  = hit.volumeId();
      nt[3]  = pos.x();
      nt[4]  = pos.y();
      nt[5]  = pos.z();
      nt[6]  = mid.x();
      nt[7]  = mid.y();
      nt[8]  = mid.z();
      nt[9]  = pca.dca();
      nt[10] = hit.time();
      nt[11] = straw.Id().getDevice();
      nt[12] = straw.Id().getSector();

      _ntup->Fill(nt);

      // Print out limited to the first few events.
      if ( ncalls < _maxFullPrint ){
	cout << "Readback hit: "
	     << evt.id().event()
	     << n++ <<  " "
	     << hit.trackId()    << "   "
	     << hit.volumeId() << " | "
	     << pca.dca()   << " "
	     << pos  << " "
	     << mom  << " "
	     << point.mag() 
	     << endl;
      }

    } // end loop over hits.
  } // end of ::analyze.


  // Count how many of this straw's nearest neighbours are hit.
  // If we have enough hits per event, it will make sense to make
  // a class to let us direct index into a list of which straws have hits.
  int ReadBack::countHitNeighbours( Straw const& straw, 
				    edm::Handle<StepPointMCCollection>& hits ){
    
    int count(0);
    vector<StrawIndex> const& nearest = straw.nearestNeighboursByIndex();
    for ( vector<int>::size_type ihit = 0;
	  ihit<nearest.size(); ++ihit ){

      StrawIndex idx = nearest[ihit];

      for( StepPointMCCollection::const_iterator 
	     i = hits->begin(),
	     e = hits->end(); i!=e ; ++i ) {
	const StepPointMC& hit = *i;
	if ( hit.strawIndex() == idx ){
	  ++count;
	  break;
	}
      }

    }
    return count;
  }
  
}


using mu2e::ReadBack;
DEFINE_FWK_MODULE(ReadBack);
