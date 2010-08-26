//
// An EDProducer Module that reads StepPointMC objects and turns them into
// StrawHit objects.
//
// $Id: MakeStrawHit_plugin.cc,v 1.3 2010/08/26 19:15:51 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/08/26 19:15:51 $
//
// Original author Rob Kutschke. Updated by Ivan Logashenko.
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

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/StrawHit.hh"
#include "ToyDP/inc/StrawHitCollection.hh"
#include "ToyDP/inc/StrawHitMCTruth.hh"
#include "ToyDP/inc/StrawHitMCTruthCollection.hh"
#include "ToyDP/inc/StrawHitMCPtr.hh"
#include "ToyDP/inc/StrawHitMCPtrCollection.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "Mu2eUtilities/inc/LinePointPCA.hh"

// Other includes.
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Vector/ThreeVector.h"

using namespace std;
using edm::Event;

namespace mu2e {

  // Utility class (structure) to hold calculated drift time for G4 hits

  class StepHit {
  public:

    int _hit_id;
    double _edep;
    double _dca;
    double _driftTime;
    double _distanceToMid;
    double _t1;
    double _t2;

    StepHit(int hit_id, double edep, double dca, double driftT, double toMid, double t1, double t2):
      _hit_id(hit_id), _edep(edep), _dca(dca), _driftTime(driftT), 
      _distanceToMid(toMid), _t1(t1), _t2(t2) { }

    // This operator is overloaded in order to time-sort the hits 
    bool StepHit::operator <(const StepHit& b) const { return (_t1 < b._t1); }

  };

  //--------------------------------------------------------------------
  //
  // 
  class MakeStrawHit : public edm::EDProducer {
  public:
    explicit MakeStrawHit(edm::ParameterSet const& pset) : 

      // Parameters
      _diagLevel(pset.getUntrackedParameter<int>("diagLevel",0)),
      _t0Sigma(pset.getUntrackedParameter<double>("t0Sigma",5.0)), // ns
      _minimumEnergy(pset.getUntrackedParameter<double>("minimumEnergy",0.0001)), // MeV
      _minimumLength(pset.getUntrackedParameter<double>("minimumLength",0.01)),   // mm
      _driftVelocity(pset.getUntrackedParameter<double>("driftVelocity",0.05)),   // mm/ns
      _driftSigma(pset.getUntrackedParameter<double>("driftSigma",0.1)),          // mm
      _minimumTimeGap(pset.getUntrackedParameter<double>("minimumTimeGap",100.0)),// ns
      _maxFullPrint(pset.getUntrackedParameter<int>("maxFullPrint",5)),
      _g4ModuleLabel(pset.getParameter<string>("g4ModuleLabel")),

      // Random number distributions
      _randGauss( createEngine( get_seed_value(pset)) ),

      _messageCategory("StrawHitMaker"){

      // Tell the framework what we make.
      produces<StrawHitCollection>();
      produces<StrawHitMCTruthCollection>();
      produces<StrawHitMCPtrCollection>();

    }
    virtual ~MakeStrawHit() { }

    virtual void beginJob(edm::EventSetup const&);
 
    void produce( edm::Event& e, edm::EventSetup const&);

  private:
    
    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Parameters
    double _t0Sigma;        // T0 spread in ns
    double _minimumEnergy;  // minimum energy deposition of G4 step 
    double _minimumLength;  // is G4Step is shorter than this, consider it a point
    double _driftVelocity;  
    double _driftSigma;
    double _minimumTimeGap; 
    string _g4ModuleLabel;  // Name of the module that made these hits.

    // Random number distributions
    CLHEP::RandGauss _randGauss;

    // A category for the error logger.
    const std::string _messageCategory;

  };

  void MakeStrawHit::beginJob(edm::EventSetup const& ){
    
  }

  void
  MakeStrawHit::produce(edm::Event& event, edm::EventSetup const&) {

    if ( _diagLevel > 0 ) cout << "MakeStrawHit: produce() begin" << endl;
      
    static int ncalls(0);
    ++ncalls;

    // Get a reference to one of the L or T trackers.
    // Throw exception if not successful.
    const Tracker& tracker = getTrackerOrThrow();

    // A container to hold the output hits.
    auto_ptr<StrawHitCollection>        strawHits(new StrawHitCollection);
    auto_ptr<StrawHitMCTruthCollection> truthHits(new StrawHitMCTruthCollection);
    auto_ptr<StrawHitMCPtrCollection>   mcptrHits(new StrawHitMCPtrCollection);

    // Ask the event to give us a handle to the requested hits.
    edm::Handle<StepPointMCCollection> points;
    event.getByLabel(_g4ModuleLabel,points);

    // Product Id of the input points.
    edm::ProductID const& id(points.id());

    // Calculate T0 for this event
    double t0 = _randGauss.fire(0.,_t0Sigma);

    // Organize hits by straws

    typedef std::map<StrawIndex,std::vector<int> > StrawHitMap;
    StrawHitMap hitmap;
    for ( int i=0; i<points->size(); ++i){
      StepPointMC const& hit = (*points)[i];
      if( hit.totalEDep()<_minimumEnergy ) continue; // Skip steps with very low energy deposition
      StrawIndex straw_id = hit.strawIndex();
      vector<int> &hits_id = hitmap[straw_id];
      hits_id.push_back(i);
    }
 
    vector<StepHit> straw_hits;

    // Loop over all straws and create StrawHits. There can be several 
    // hits per straw if they are separated by time. The general algorithm 
    // is as follows: calculate signal time for each G4step, order them 
    // in time and look for gaps. If gap exceeds _minimumTimeGap create 
    // separate hit.

    for(StrawHitMap::const_iterator istraw = hitmap.begin(); istraw != hitmap.end(); ++istraw) {

      if ( ncalls < _maxFullPrint && _diagLevel > 2 ) {
        cout << "MakeStrawHit: straw ID=" << istraw->first 
             << ": number of G4 step hits " << istraw->second.size() 
             << endl;
      }

      // Get the straw information, also by reference.
      StrawIndex straw_id = istraw->first;
      Straw const&      straw = tracker.getStraw(straw_id);
      CLHEP::Hep3Vector const& mid   = straw.getMidPoint();
      CLHEP::Hep3Vector const& w     = straw.getDirection();
      double strawHalfLength         = straw.getHalfLength();

      // Prepare info for hit creation
      straw_hits.clear();

      // Loop over all hits found for this straw

      vector<int> const& ihits = istraw->second;

      for( int i=0; i<ihits.size(); i++ ) {

        int hitRef = ihits[i];
        StepPointMC const& hit = (*points)[hitRef];
        CLHEP::Hep3Vector  const& pos = hit.position();
        CLHEP::Hep3Vector  const& mom = hit.momentum();
        double length  = hit.stepLength();
        double edep    = hit.totalEDep();
        double hitTime = hit.time();

        if ( ncalls < _maxFullPrint && _diagLevel > 2 ) {
          cout << "MakeStrawHit: Hit #" << i << " : length=" << length 
               << " energy=" << edep << " time=" << hitTime 
               << endl;
        }
        
        // Calculate the drift distance from this step. 
        double hit_dca;
        CLHEP::Hep3Vector hit_pca;

        if( length < _minimumLength ) {

          // If step length is very small, consider it a point

          LinePointPCA pca(mid, w, pos);
          hit_dca = pca.dca();
          hit_pca = pca.pca();
          
        } else {

          // Step is not a point. Calculate the distance between two lines. 

          TwoLinePCA pca( mid, w, pos, mom);
          CLHEP::Hep3Vector const& p2 = pca.point2();

          if( (pos-p2).mag()<=length && (pos-p2).dot(mom)<=0 ) {

            // If the point of closest approach is within the step and wire - thats it.
            hit_dca = pca.dca();
            hit_pca = pca.point1();

          } else {
            
            // The point of closest approach is not within the step. In this case
            // the closes distance should be calculated from the ends

            LinePointPCA pca1(mid, w, pos);
            LinePointPCA pca2(mid, w, pos+mom.unit()*length);
            if( pca1.dca() < pca2.dca() ) {
              hit_dca = pca1.dca();
              hit_pca = pca1.pca();
            } else {
              hit_dca = pca2.dca();
              hit_pca = pca2.pca();
            } 
          
          }          

        } // drift distance calculation

        // Calculate signal time. It is Geant4 time + signal propagation time
        // t1 is signal time at positive end (along w vector), 
        // t2 - at negative end (opposite to w vector)

        const double signalVelocity = 299.792458; // mm/ns

        double driftTime = (hit_dca + _randGauss.fire(0.,_driftSigma))/_driftVelocity;
        double distanceToMiddle = (hit_pca-mid).dot(w);
        double hit_t1 = t0 + hitTime + driftTime + (strawHalfLength-distanceToMiddle)/signalVelocity;
        double hit_t2 = t0 + hitTime + driftTime + (strawHalfLength+distanceToMiddle)/signalVelocity;

        straw_hits.push_back(StepHit(hitRef,edep,hit_dca,driftTime,distanceToMiddle,hit_t1,hit_t2));

      } // loop over hits

      // Now that we calculated estimated signal time for all G4Steps, we can analyze
      // the time structure and create StrawHits

      // First we need to sort StepHits according to t1 time
      sort(straw_hits.begin(), straw_hits.end() );

      if ( ncalls < _maxFullPrint && _diagLevel > 2 ) {
        for( int i=0; i<straw_hits.size(); i++ ) {
          cout << "MakeStrawHit: StepHit #" << straw_hits[i]._hit_id 
               << " DCA=" << straw_hits[i]._dca
               << " driftT=" << straw_hits[i]._driftTime
               << " distToMid=" << straw_hits[i]._distanceToMid
               << " t1=" << straw_hits[i]._t1
               << " t2=" << straw_hits[i]._t2
               << " edep=" << straw_hits[i]._edep
               << " t0=" << t0
               << endl;
        }
      }

      // Now loop over all StepHits and create StrawHits as needed

      if( straw_hits.size()<1 ) continue; // This should never be needed. Added for safety.

      double digi_time   = straw_hits[0]._t1;
      double digi_dt     = straw_hits[0]._t2 - straw_hits[0]._t1;
      double digi_edep   = straw_hits[0]._edep;
      double digi_driftT = straw_hits[0]._driftTime;
      double digi_toMid  = straw_hits[0]._distanceToMid;
      double digi_dca    = straw_hits[0]._dca;
      StrawHitMCPtr mcptr;
      mcptr.push_back(DPIndex(id,straw_hits[0]._hit_id));

      for( int i=1; i<straw_hits.size(); i++ ) {
        if( (straw_hits[i]._t1-straw_hits[i-1]._t1) > _minimumTimeGap ) {
          // The is bit time gap - save current data as a hit...
          strawHits->push_back(StrawHit(straw_id,digi_time,digi_dt,digi_edep));
          truthHits->push_back(StrawHitMCTruth(t0,digi_driftT,digi_dca,digi_toMid));
          mcptrHits->push_back(mcptr);
          // ...and create new hit
          mcptr.clear();
          mcptr.push_back(DPIndex(id,straw_hits[i]._hit_id));
          digi_time   = straw_hits[i]._t1;
          digi_dt     = straw_hits[i]._t2 - straw_hits[i]._t1;
          digi_edep   = straw_hits[i]._edep;
          digi_driftT = straw_hits[i]._driftTime;
          digi_toMid  = straw_hits[i]._distanceToMid;
          digi_dca    = straw_hits[i]._dca;
        } else {
          // Append existing hit
          digi_edep += straw_hits[i]._edep;
          mcptr.push_back(DPIndex(id,straw_hits[i]._hit_id));
        }
      }

      strawHits->push_back(StrawHit(straw_id,digi_time,digi_dt,digi_edep));
      truthHits->push_back(StrawHitMCTruth(t0,digi_driftT,digi_dca,digi_toMid));
      mcptrHits->push_back(mcptr);

    }
    
    if ( ncalls < _maxFullPrint && _diagLevel > 2 ) {
      cout << "MakeStrawHit: Total number of hit straws = " << strawHits->size() << endl;
      for( int i=0; i<strawHits->size(); ++i ) {
        cout << "MakeStrawHit: Straw #" << (*strawHits)[i].strawIndex() 
             << " time="  << (*strawHits)[i].time()
             << " dt="    << (*strawHits)[i].dt()
             << " edep="  << (*strawHits)[i].energyDep()
             << " nsteps="<< (*mcptrHits)[i].size()
             << endl;
      }
    }

    // Add the output hit collection to the event
    event.put(strawHits);
    event.put(truthHits);
    event.put(mcptrHits);

    if ( _diagLevel > 0 ) cout << "MakeStrawHit: produce() end" << endl;

  } // end of ::analyze.
  
}

using mu2e::MakeStrawHit;
DEFINE_FWK_MODULE(MakeStrawHit);
