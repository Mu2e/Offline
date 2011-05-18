//
// An EDProducer Module that reads StepPointMC objects and turns them into
// StrawHit objects.
//
// $Id: MakeDriftCellHit_module.cc,v 1.4 2011/05/18 15:47:40 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/18 15:47:40 $
//
// Original author G.F. Tassielli. Class derived by MakeStrawHit
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

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "ToyDP/inc/StrawHit.hh"
#include "ToyDP/inc/StrawHitCollection.hh"
#include "ToyDP/inc/StrawHitMCTruth.hh"
#include "ToyDP/inc/StrawHitMCTruthCollection.hh"
#include "ToyDP/inc/DPIndexVector.hh"
#include "ToyDP/inc/DPIndexVectorCollection.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "Mu2eUtilities/inc/LinePointPCA.hh"

// Other includes.
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Vector/ThreeVector.h"

using namespace std;
using art::Event;

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

    StepHit(int hit_id, double edep, double dca, double driftT, double toMid, double t1):
      _hit_id(hit_id), _edep(edep), _dca(dca), _driftTime(driftT),
      _distanceToMid(toMid), _t1(t1) { }

    // This operator is overloaded in order to time-sort the hits
    bool operator <(const StepHit& b) const { return (_t1 < b._t1); }

  };

  //--------------------------------------------------------------------
  //
  //
  class MakeDriftCellHit : public art::EDProducer {
  public:
    explicit MakeDriftCellHit(fhicl::ParameterSet const& pset) :

      // Parameters
      _diagLevel(pset.get<int>("diagLevel",0)),
      _maxFullPrint(pset.get<int>("maxFullPrint",5)),
      _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
      _t0Sigma(pset.get<double>("t0Sigma",5.0)), // ns
      _timetodist(pset.get<double>("timetodist",0.0/*149.8962*/)), // mm/ns
      _distSigma(pset.get<double>("distSigma",0.0/*80.*/)), // mm
      _minimumEnergy(pset.get<double>("minimumEnergy",0.0001)), // MeV
      _minimumLength(pset.get<double>("minimumLength",0.005/*0.01*/)),   // mm
      _driftVelocity(pset.get<double>("driftVelocity",0.035/*0.05*/)),   // mm/ns
      _driftSigma(pset.get<double>("driftSigma",0.2/*0.1*/)),          // mm
      _minimumTimeGap(pset.get<double>("minimumTimeGap",500.0/*100.0*/)),// ns depends by the cell dimension it will change after
      _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),

      // Random number distributions
      _gaussian( createEngine( get_seed_value(pset)) ),

      _messageCategory("DriftCellHitMaker"/*"StrawHitMaker"*/){

      // Tell the framework what we make.
      produces<StrawHitCollection>();
      produces<StrawHitMCTruthCollection>();
      produces<DPIndexVectorCollection>("StrawHitMCPtr");

    }
    virtual ~MakeDriftCellHit() { }

    virtual void beginJob();

    void produce( art::Event& e);

  private:

    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;

    // Parameters
    double _t0Sigma;        // T0 spread in ns
    double _timetodist;     // const to convert delata t in delat z along the wire in mm/ns
    double _distSigma;      // sigma of dealta z in mm
    double _minimumEnergy;  // minimum energy deposition of G4 step
    double _minimumLength;  // is G4Step is shorter than this, consider it a point
    double _driftVelocity;
    double _driftSigma;
    double _minimumTimeGap;
    string _g4ModuleLabel;  // Name of the module that made these hits.

    // Random number distributions
    CLHEP::RandGaussQ _gaussian;

    // A category for the error logger.
    const std::string _messageCategory;

  };

  void MakeDriftCellHit::beginJob(){

  }

  void
  MakeDriftCellHit::produce(art::Event& event) {

    if ( _diagLevel > 0 ) cout << "MakeDriftCellHit: produce() begin" << endl;

    static int ncalls(0);
    ++ncalls;

    // Get a reference to one of the I trackers.
    // Throw exception if not successful.
    //const Tracker& tracker = getTrackerOrThrow();
    const Tracker& tracker = getTrackerOrThrow();
    //const ITracker& tracker     = static_cast<const ITracker&>(tracker);
    //GeomHandle<ITracker> tracker;

    // A container to hold the output hits.
    auto_ptr<StrawHitCollection>        strawHits(new StrawHitCollection);
    auto_ptr<StrawHitMCTruthCollection> truthHits(new StrawHitMCTruthCollection);
    auto_ptr<DPIndexVectorCollection>   mcptrHits(new DPIndexVectorCollection);

    // Ask the event to give us a handle to the requested hits.
    art::Handle<StepPointMCCollection> points;
    event.getByLabel(_g4ModuleLabel,_trackerStepPoints,points);

    // Product Id of the input points.
    art::ProductID const& id(points.id());

    // Calculate T0 for this event
    double t0 = _gaussian.fire(0.,_t0Sigma);

    // Organize hits by cells

    //typedef std::map<StepPointMC::VolumeId_type,std::vector<int> > DriftCellHitMap;
    typedef std::map<StrawIndex,std::vector<int> > DriftCellHitMap;
    DriftCellHitMap hitmap;
    for ( size_t i=0; i<points->size(); ++i){
      StepPointMC const& hit = (*points)[i];
      if( hit.totalEDep()<_minimumEnergy ) continue; // Skip steps with very low energy deposition
      //StepPointMC::VolumeId_type dcell_id = hit.volumeId();
      StrawIndex dcell_id = hit.strawIndex();
      vector<int> &hits_id = hitmap[dcell_id];
      hits_id.push_back(i);
    }

    vector<StepHit> dcell_hits;

    // Loop over all drift cells and create StrawHits. There can be several
    // hits per cell if they are separated by time. The general algorithm
    // is as follows: calculate signal time for each G4step, order them
    // in time and look for gaps. If gap exceeds _minimumTimeGap create
    // separate hit.

    for(DriftCellHitMap::const_iterator idcell = hitmap.begin(); idcell != hitmap.end(); ++idcell) {

      if ( ncalls < _maxFullPrint && _diagLevel > 2 ) {
        cout << "MakeDriftCellHit: cell ID=" << idcell->first
             << ": number of G4 step hits " << idcell->second.size()
             << endl;
      }

      // Get the straw information, also by reference.
//      StepPointMC::VolumeId_type dcell_id = idcell->first;
//      tracker->getCellGeometryHandle()->SelectCellDet(dcell_id);
//      boost::shared_ptr<Cell>   cell = tracker->getCellGeometryHandle()->GetITCell();
//
//      CLHEP::Hep3Vector const& mid   = cell->getMidPoint();
//      CLHEP::Hep3Vector const& w     = cell->getDirection();
//      double strawHalfLength         = cell->getHalfLength();

      StrawIndex dcell_id = idcell->first;
      //const Straw&  cell = tracker.getStraw(dcell_id);
      const Cell&    cell = static_cast<const Cell&>( tracker.getStraw(dcell_id) );
      CLHEP::Hep3Vector const& mid   = cell.getMidPoint();
      CLHEP::Hep3Vector const& w     = cell.getDirection();
      double strawHalfLength         = cell.getHalfLength();

//      StrawIndex dcell_id = idcell->first;
//      itracker->getCellGeometryHandle()->SelectCellDet(dcell_id.asUint());
//      boost::shared_ptr<Cell>   cell = itracker->getCellGeometryHandle()->GetITCell();
//
//      CLHEP::Hep3Vector const& mid   = cell->getMidPoint();
//      CLHEP::Hep3Vector const& w     = cell->getDirection();
//      double strawHalfLength         = cell->getHalfLength();


      //_minimumTimeGap                = 0.5*(cell->getDetail()->InscribedCircleRadius()+cell->getDetail()->CirumscribedRadius())/_driftVelocity;
      _minimumTimeGap                = 0.5*(cell.getDetail()->InscribedCircleRadius()+cell.getDetail()->CirumscribedRadius())/_driftVelocity;

      // Prepare info for hit creation
      dcell_hits.clear();

      // Loop over all hits found for this straw

      vector<int> const& ihits = idcell->second;

      for( size_t i=0; i<ihits.size(); i++ ) {

        int hitRef = ihits[i];
        StepPointMC const& hit = (*points)[hitRef];
        CLHEP::Hep3Vector  const& pos = hit.position();
        CLHEP::Hep3Vector  const& mom = hit.momentum();
        double length  = hit.stepLength();
        double edep    = hit.totalEDep();
        double hitTime = hit.time();

        if ( ncalls < _maxFullPrint && _diagLevel > 2 ) {
          cout << "MakeDriftCellHit: Hit #" << i << " : length=" << length
               << " energy=" << edep << " time=" << hitTime
               << endl;
        }

        // Calculate the drift distance from this step.
        double hit_dca,sign;
        CLHEP::Hep3Vector hit_pca;

        if( length < _minimumLength ) {

          // If step length is very small, consider it a point

          LinePointPCA pca(mid, w, pos);
          hit_dca = pca.dca();
          hit_pca = pca.pca();
	  continue;
        } else {

          // Step is not a point. Calculate the distance between two lines.

          TwoLinePCA pca( mid, w, pos, mom);
          CLHEP::Hep3Vector const& p2 = pca.point2();

	  //          if( (pos-p2).mag()<=length && (pos-p2).dot(mom)<=0 ) {

            // If the point of closest approach is within the step and wire - thats it.
            hit_dca = pca.dca();
            hit_pca = pca.point1();
	    sign=w.cross(mom).dot(hit_pca-p2);
	    /*
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
	    */
        } // drift distance calculation

        // Calculate signal time. It is Geant4 time + signal propagation time
        // t1 is signal time at positive end (along w vector),
        // (Not used for ITracker) t2 - at negative end (opposite to w vector)

        const double signalVelocity = 299.792458; // mm/ns

        double driftTime = (hit_dca + _gaussian.fire(0.,_driftSigma))/_driftVelocity;
        double distanceToMiddle = (hit_pca-mid).dot(w);
        double hit_t1 = t0 + hitTime + driftTime + (strawHalfLength-distanceToMiddle)/signalVelocity;
        //double hit_t2 = -9999.9;//t0 + hitTime + driftTime + (strawHalfLength+distanceToMiddle)/signalVelocity;

        dcell_hits.push_back(StepHit(hitRef,edep,hit_dca*(sign>0?1.:-1.),driftTime,distanceToMiddle,hit_t1));

      } // loop over hits

      // Now that we calculated estimated signal time for all G4Steps, we can analyze
      // the time structure and create StrawHits

      // First we need to sort StepHits according to t1 time
      sort(dcell_hits.begin(), dcell_hits.end() );

      if ( ncalls < _maxFullPrint && _diagLevel > 2 ) {
        for( size_t i=0; i<dcell_hits.size(); i++ ) {
          cout << "MakeDriftCellHit: StepHit #" << dcell_hits[i]._hit_id
               << " DCA=" << dcell_hits[i]._dca
               << " driftT=" << dcell_hits[i]._driftTime
               << " distToMid=" << dcell_hits[i]._distanceToMid
               << " t1=" << dcell_hits[i]._t1
               /*<< " t2=" << dcell_hits[i]._t2*/
               << " edep=" << dcell_hits[i]._edep
               << " t0=" << t0
               << endl;
        }
      }

      // Now loop over all StepHits and create StrawHits as needed

      if( dcell_hits.size()<1 ) continue; // This should never be needed. Added for safety.

      double digi_time   = dcell_hits[0]._t1;
      double digi_t2     = 0.0;//dcell_hits[0]._t2;
      double digi_edep   = dcell_hits[0]._edep;
      double digi_driftT = dcell_hits[0]._driftTime;
      double digi_toMid  = dcell_hits[0]._distanceToMid;
      double digi_dca    = dcell_hits[0]._dca;
      double deltadigitime;
      DPIndexVector mcptr;
      mcptr.push_back(DPIndex(id,dcell_hits[0]._hit_id));

      for( size_t i=1; i<dcell_hits.size(); i++ ) {
        if( (dcell_hits[i]._t1-dcell_hits[i-1]._t1) > _minimumTimeGap ) {
          // The is bit time gap - save current data as a hit...
          strawHits->push_back(StrawHit(/*StrawIndex(*/dcell_id/*)*/,digi_time,digi_t2/*digi_t2-digi_time*/,digi_edep));
          truthHits->push_back(StrawHitMCTruth(t0,digi_driftT,digi_dca,digi_toMid));
          mcptrHits->push_back(mcptr);
          // ...and create new hit
          mcptr.clear();
          mcptr.push_back(DPIndex(id,dcell_hits[i]._hit_id));
          digi_time   = dcell_hits[i]._t1;
          //digi_t2     = dcell_hits[i]._t2;
          digi_edep   = dcell_hits[i]._edep;
          digi_driftT = dcell_hits[i]._driftTime;
          digi_toMid  = dcell_hits[i]._distanceToMid;
          digi_dca    = dcell_hits[i]._dca;
        } else {
          // Append existing hit
	  //if( digi_t2 > dcell_hits[i]._t2 ) digi_t2 = dcell_hits[i]._t2;
          digi_edep += dcell_hits[i]._edep;
          mcptr.push_back(DPIndex(id,dcell_hits[i]._hit_id));
        }
      }
      //deltadigitime=(digi_t2-digi_time)+_gaussian.fire(0.,_distSigma/_timetodist);
      deltadigitime=0.0;
      strawHits->push_back(StrawHit(/*StrawIndex(*/dcell_id/*)*/,digi_time,deltadigitime/*deltadigitime*/,digi_edep));
      truthHits->push_back(StrawHitMCTruth(t0,digi_driftT,digi_dca,digi_toMid));
      mcptrHits->push_back(mcptr);

    }

    if ( ncalls < _maxFullPrint && _diagLevel > 2 ) {
      cout << "MakeDriftCellHit: Total number of hit cells = " << strawHits->size() << endl;
      for( size_t i=0; i<strawHits->size(); ++i ) {
        cout << "MakeDriftCellHit: Cell #" << (*strawHits)[i].strawIndex()
             << " time="  << (*strawHits)[i].time()
             /*<< " dt="    << (*strawHits)[i].dt()*/
             << " edep="  << (*strawHits)[i].energyDep()
             << " nsteps="<< (*mcptrHits)[i].size()
             << endl;
      }
    }

    // Add the output hit collection to the event
    event.put(strawHits);
    event.put(truthHits);
    event.put(mcptrHits,"StrawHitMCPtr");

    if ( _diagLevel > 0 ) cout << "MakeDriftCellHit: produce() end" << endl;

  } // end of ::analyze.

}

using mu2e::MakeDriftCellHit;
DEFINE_ART_MODULE(MakeDriftCellHit);
