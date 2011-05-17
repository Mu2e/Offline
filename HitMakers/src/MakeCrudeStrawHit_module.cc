//
// An EDProducer Module that reads StepPointMC objects and turns them into
// CrudeStrawHit objects.
//
// $Id: MakeCrudeStrawHit_module.cc,v 1.1 2011/05/17 16:30:14 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 16:30:14 $
//
// Original author Rob Kutschke
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
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
#include "TH2F.h"
#include "TNtuple.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "LTrackerGeom/inc/CrudeStrawHitCollection.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"

// Other includes.
#include "CLHEP/Random/RandGaussQ.h"

using namespace std;

namespace mu2e {

  //--------------------------------------------------------------------
  //
  // 
  class MakeCrudeStrawHit : public art::EDProducer {
  public:
    explicit MakeCrudeStrawHit(fhicl::ParameterSet const& pset) : 

      // Parameters
      _diagLevel(pset.get<int>("diagLevel",0)),
      _maxFullPrint(pset.get<int>("maxFullPrint",5)),

      // Random number distributions
      _gaussian( createEngine( get_seed_value(pset)) ),

      // Histograms
      _hTime(0),
      _hDriftDist(0),
      _hCheckPointRadius(0),

      _messageCategory("HitMaker"){

      // A place holder.
      produces<CrudeStrawHitPData>();

    }
    virtual ~MakeCrudeStrawHit() { }

    virtual void beginJob(art::EventSetup const&);
 
    void produce( art::Event& e, art::EventSetup const&);

  private:
    
    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    CLHEP::RandGaussQ _gaussian;

    // Pointers to diagnostic histograms.
    TH1F* _hTime;
    TH1F* _hDriftDist;
    TH1F* _hCheckPointRadius;

    // A category for the error logger.
    const std::string _messageCategory;

  };

  void MakeCrudeStrawHit::beginJob(art::EventSetup const& ){

    // Create histograms if diagnostics are enabled.
    if ( _diagLevel > 0 ){

      art::ServiceHandle<art::TFileService> tfs;

      _hTime      = tfs->make<TH1F>( "hTime", "Pulse Height;(ns)",              100,  0.,  2000. );
      _hDriftDist = tfs->make<TH1F>( "hDriftDist", "Crude Drift Distance;(mm)", 100,  0.,     3.  );
      _hCheckPointRadius = tfs->make<TH1F>( "hCheckPointRadius",  "Radius of Reference point; (mm)",
                                            100, 2.25, 2.75 );
    }

  }

  void
  MakeCrudeStrawHit::produce(art::Event& event, art::EventSetup const&) {

    static int ncalls(0);
    ++ncalls;

    // Geometry for the LTracker.
    GeomHandle<LTracker> ltracker;

    // A container to hold the output hits.
    auto_ptr<CrudeStrawHitPData> crudeHits(new CrudeStrawHitPData);

    // Instance name of the module that created the hits of interest;
    static const string creatorName("g4run");
    static const string collectionName("tracker");

    // Ask the event to give us a handle to the requested hits.
    art::Handle<StepPointMCCollection> points;
    event.getByLabel(creatorName,collectionName,points);

    // Product Id of the input points.
    art::ProductID const& id(points.id());

    for ( size_t i=0; i<points->size(); ++i){
      
      // Aliases (references), used for readability.
      StepPointMC const& hit = (*points)[i];
      CLHEP::Hep3Vector  const& pos = hit.position();
      CLHEP::Hep3Vector  const& mom = hit.momentum();
      
      // Get the straw information, also by reference.
      Straw const&      straw = ltracker->getStraw(hit.strawIndex());
      CLHEP::Hep3Vector const& mid   = straw.getMidPoint();
      CLHEP::Hep3Vector const& w     = straw.getDirection();

      // Compute straight line approximation of the drift distance.
      TwoLinePCA pca( mid, w, pos, mom);
      double dcaTrue = pca.dca();

      // Move these somewhere smart.
      const static double driftVelocity = 0.05;  // mm/ns
      const static double sigma = 0.100;         // mm

      // Compute drift time and smeared drift distance.
      double time = dcaTrue/driftVelocity;
      double dca  = dcaTrue + _gaussian.fire(0.,sigma);

      // Add to the output collection.
      crudeHits->push_back( CrudeStrawHit( straw.Index(), 
                                           dca, 
                                           time, 
                                           sigma, 
                                           hit.eDep(), 
                                           CrudeStrawHit::stepPointMC,
                                           DPIndex(id,i),
                                           dcaTrue,
                                           &event)
                            );

      // Fill diagnostic histograms.
      if ( _diagLevel > 0 ) {

        // Check the radius of the reference point in the local
        // coordinates of the straw.  It should be 2.5 mm.
        double s = w.dot(pos-mid);
        CLHEP::Hep3Vector point = pos - (mid + s*w);

        // I don't understand the distribution of the time variable.
        // I want it to be the time from the start of the spill.
        // It appears to be the time since start of tracking.

        _hTime->Fill(hit.time());
        _hCheckPointRadius->Fill(point.mag());
        _hDriftDist->Fill(pca.dca());
      }


    } // end loop over StepPointMC's

    // Diagnostic printout.
    if ( ncalls < _maxFullPrint && _diagLevel > 2){
      for ( unsigned int i=0; i<crudeHits->size(); ++i){
        cout << event.id().event() << " | "   << (*crudeHits)[i] << endl;
      }
    }

    // All done.  Add the output to the event.
    event.put(crudeHits);

    
  } // end of ::analyze.
  
}


using mu2e::MakeCrudeStrawHit;
DEFINE_ART_MODULE(MakeCrudeStrawHit);
