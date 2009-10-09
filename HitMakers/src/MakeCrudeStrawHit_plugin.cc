//
// An EDProducer Module that reads StepPointMC objects and turns them into
// crude straw hit objects.
//
// $Id: MakeCrudeStrawHit_plugin.cc,v 1.1 2009/10/09 13:31:32 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/10/09 13:31:32 $
//
// Original author Rob Kutschke
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

// Framework includes.
#include "FWCore/Framework/interface/EDProducer.h"
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
#include "TH2F.h"
#include "TNtuple.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/CrudeStrawHitCollection.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"

// Other includes.
#include "CLHEP/Random/RandGauss.h"


using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e {

  //--------------------------------------------------------------------
  //
  // 
  class MakeCrudeStrawHit : public edm::EDProducer {
  public:
    explicit MakeCrudeStrawHit(edm::ParameterSet const& pset) : 
      _diagLevel(0),
      _maxFullPrint(pset.getUntrackedParameter<int>("maxFullPrint",5)),
      _hRadius(0),
      _hTime(0),
      _hMultiplicity(0),
      _hDriftDist(0),
      _messageCategory("HitMaker"){

      // A place holder.
      produces<CrudeStrawHitCollection>();

    }
    virtual ~MakeCrudeStrawHit() { }

    virtual void beginJob(edm::EventSetup const&);
 
    void produce( edm::Event& e, edm::EventSetup const&);

  private:
    

    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

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

  void MakeCrudeStrawHit::beginJob(edm::EventSetup const& ){

    // Get access to the TFile service.
    edm::Service<edm::TFileService> tfs;

    // Create some 1D histograms.
    _hRadius       = tfs->make<TH1F>( "hRadius", "Radius of Hits;(mm)",          100,  0., 1000. );
    _hTime         = tfs->make<TH1F>( "hTime", "Pulse Height;(ns)",              100,  0.,  100. );
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

  void
  MakeCrudeStrawHit::produce(edm::Event& evt, edm::EventSetup const&) {

    static int ncalls(0);
    ++ncalls;

    // Master geometry for the LTracker.
    GeomHandle<LTracker> ltracker;

    // Hold the output hits.
    auto_ptr<CrudeStrawHitCollection> crudeHits(new CrudeStrawHitCollection);

    // Instance name of the module that created the hits of interest;
    static const string creatorName("g4run");

    // Ask the event to give us a "handle" to the requested hits.
    edm::Handle<StepPointMCCollection> hits;
    evt.getByLabel(creatorName,hits);
    
    // Fill histogram with number of hits per event.
    _hMultiplicity->Fill(hits->size());

    // ntuple buffer.
    float nt[13];


    for ( int i=0; i<hits->size(); ++i){
      
      // Aliases, used for readability.
      StepPointMC const& hit = (*hits)[i];
      Hep3Vector  const& pos = hit.position();
      Hep3Vector  const& mom = hit.momentum();
      
      // Get the straw information.
      Straw const&      straw = ltracker->getStraw( hit.volumeId() );
      Hep3Vector const& mid   = straw.getMidPoint();
      Hep3Vector const& w     = straw.getDirection();

      // Count how many nearest neighbours are also hit.
      int nNeighbours = countHitNeighbours( straw, hits );

      // Compute an estimate of the drift distance.
      TwoLinePCA pca( mid, w, pos, mom);
      double dcaTrue = pca.dca();

      // Move these somewhere smart.
      const static double driftVelocity = 0.05;  // mm/ns
      const static double sigma = 0.100;         // mm

      // Compute drift time and smeared drift distance.
      double time = dcaTrue/driftVelocity;
      double dca  = dcaTrue + RandGauss::shoot(0.,sigma);

      crudeHits->push_back( CrudeStrawHit( straw.Index(), 
					   dca, 
					   time, 
					   sigma, 
					   hit.eDep(), 
					   i, 
					   dcaTrue)
			    );

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
      if ( ncalls < _maxFullPrint && _diagLevel > 0){
	cout << "Readback hit: "
	     << evt.id().event() << " "
	     << i                <<  " "
	     << hit.trackId()  << "   "
	     << hit.volumeId() << " | "
	     << pca.dca()      << " "
	     << pos  << " "
	     << mom  << " "
	     << point.mag()  << " " 
	     << point.z()
	     << endl;

      }

    } // end loop over hits.

    for ( unsigned int i=0; i<crudeHits->size(); ++i){
      cout << evt.id().event() << " | "   << (*crudeHits)[i] << endl;
    }

    evt.put(crudeHits);
    
  } // end of ::analyze.


  // Count how many of this straw's nearest neighbours are hit.
  // If we have enough hits per event, it will make sense to make
  // a class to let us direct index into a list of which straws have hits.
  int MakeCrudeStrawHit::countHitNeighbours( Straw const& straw, 
				    edm::Handle<StepPointMCCollection>& hits ){
    
    int count(0);
    vector<int> const& nearest = straw.nearestNeighboursByIndex();
    for ( vector<int>::size_type ihit =0;
	  ihit<nearest.size(); ++ihit ){

      int idx = nearest[ihit];

      for( StepPointMCCollection::const_iterator 
	     i = hits->begin(),
	     e = hits->end(); i!=e ; ++i ) {
	const StepPointMC& hit = *i;
	if ( hit.volumeId() == idx ){
	  ++count;
	  break;
	}
      }

    }
    return count;
  }
  
}


using mu2e::MakeCrudeStrawHit;
DEFINE_FWK_MODULE(MakeCrudeStrawHit);
