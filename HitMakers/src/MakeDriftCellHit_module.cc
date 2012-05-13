//
// An EDProducer Module that reads StepPointMC objects and turns them into
// StrawHit objects.
//
// $Id: MakeDriftCellHit_module.cc,v 1.17 2012/05/13 21:23:15 ignatov Exp $
// $Author: ignatov $
// $Date: 2012/05/13 21:23:15 $
//
// Original author G.F. Tassielli. Class derived by MakeStrawHit
//

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ITrackerGeom/inc/ITracker.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "Mu2eUtilities/inc/LinePointPCA.hh"
#include "SeedService/inc/SeedService.hh"

// Includes from art and its tool chain.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

// Other includes.
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Vector/ThreeVector.h"

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

using namespace std;
//using art::Event;

namespace mu2e {

  // Utility class (structure) to hold calculated drift time for G4 hits

  class StepHit {
    typedef art::Handle<StepPointMCCollection> const* PHandle;

  public:

    art::Ptr<StepPointMC> _ptr;
    double _edep;
    double _dca;
    double _driftTime;
    double _distanceToMid;
    double _t1;

    StepHit( art::Ptr<StepPointMC> const& ptr, double edep, double dca, double driftT, double toMid, double t1):
      _ptr(ptr), _edep(edep), _dca(dca), _driftTime(driftT),
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
      _timetodist(pset.get<double>("timetodist",0.0/*149.8962*/)), // mm/ns
      _distSigma(pset.get<double>("distSigma",0.0/*80.*/)), // mm
      _minimumEnergy(pset.get<double>("minimumEnergy",0.00001)), // MeV
      _minimumLength(pset.get<double>("minimumLength",0.005/*0.01*/)),   // mm
      _minimumTimeGap(pset.get<double>("minimumTimeGap",500.0/*100.0*/)),// ns depends by the cell dimension it will change after
      _nzHittingZone(pset.get<int>("nzHittingZone",1)),
      _zZoneLimits(pset.get< vector<double> >("zZoneLimits")),
      _zZoneActive(pset.get< vector<bool> >("zZoneActive")),
      _zZoneMinLiveTime(pset.get< vector<double> >("zZoneMinLiveTime")),
      _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),

      // Random number distributions
      _gaussian( createEngine( art::ServiceHandle<SeedService>()->getSeed() ) ),

      _messageCategory("DriftCellHitMaker"/*"StrawHitMaker"*/),

      // Control some information messages.
      _firstEvent(true){

      // Tell the framework what we make.
      produces<StrawHitCollection>();
      produces<StrawHitMCTruthCollection>();
      produces<PtrStepPointMCVectorCollection>("StrawHitMCPtr");

      if ( _nzHittingZone>1 ) _zHittingZonePresent=true;
      else _zHittingZonePresent=false;

    }
    virtual ~MakeDriftCellHit() { }

    virtual void beginJob();

    void produce( art::Event& e);

  private:

    typedef std::map<StrawIndex,std::vector<art::Ptr<StepPointMC> > > DriftCellHitMap;

    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;

    // Parameters
    double _timetodist;     // const to convert delata t in delat z along the wire in mm/ns
    double _distSigma;      // sigma of dealta z in mm
    double _minimumEnergy;  // minimum energy deposition of G4 step
    double _minimumLength;  // is G4Step is shorter than this, consider it a point
    double _driftVelocity;
    double _driftSigma;
    double _minimumTimeGap;
    int    _nzHittingZone;
    bool   _zHittingZonePresent;
    std::vector<double> _zZoneLimits;
    std::vector<bool>   _zZoneActive;
    std::vector<double> _zZoneMinLiveTime;
    string _g4ModuleLabel;  // Name of the module that made these hits.

    // Random number distributions
    CLHEP::RandGaussQ _gaussian;

    // A category for the error logger.
    const std::string _messageCategory;

    // Give some informationation messages only on the first event.
    bool _firstEvent;

    void fillHitMap ( art::Event const& event, DriftCellHitMap& hitmap );

  };

  void MakeDriftCellHit::beginJob(){

  }


  // Find StepPointMCs in the event and use them to fill the hit map.
  void MakeDriftCellHit::fillHitMap ( art::Event const& event, DriftCellHitMap& hitmap){

    // This selector will select only data products with the given instance name.
    art::ProductInstanceNameSelector selector("tracker");

    typedef std::vector< art::Handle<StepPointMCCollection> > HandleVector;

    // Get all of the tracker StepPointMC collections from the event:
    HandleVector stepsHandles;
    event.getMany( selector, stepsHandles);

    // Informational message on the first event.
    if ( _firstEvent ) {
      _firstEvent = false;
      mf::LogInfo log(_messageCategory);
      log << "MakeDriftCellHit::fillHitMap will uses StepPointMCs from: \n";
      for ( HandleVector::const_iterator i=stepsHandles.begin(), e=stepsHandles.end();
            i != e; ++i ){

        art::Provenance const& prov(*(i->provenance()));
        log  << "   " << prov.branchName() << "\n";

      }
    }

    // Populate hitmap from stepsHandles.
    for ( HandleVector::const_iterator i=stepsHandles.begin(), e=stepsHandles.end();
          i != e; ++i ){

      art::Handle<StepPointMCCollection> const& handle(*i);
      StepPointMCCollection const& steps(*handle);

      int index(0);
      for ( StepPointMCCollection::const_iterator j = steps.begin(), je=steps.end();
            j != je; ++j, ++index){

        StepPointMC const& step(*j);
        if( step.totalEDep()<_minimumEnergy ) continue; // Skip steps with very low energy deposition
        StrawIndex straw_id = step.strawIndex();

        // The main event for this method.
        hitmap[straw_id].push_back( art::Ptr<StepPointMC>(handle,index) );
      }
    }

  }   // end MakeStrawHit::fillHitMap


  void
  MakeDriftCellHit::produce(art::Event& event) {
    //cout<<"Event "<<event.event()<<endl;
    if ( _diagLevel > 0 ) cout << "MakeDriftCellHit: produce() begin" << endl;

    static int ncalls(0);
    ++ncalls;

    // Get a reference to one of the I trackers.
    // Throw exception if not successful.
    //const Tracker& tracker = getTrackerOrThrow();
    const Tracker& tracker = getTrackerOrThrow();
    //const ITracker& tracker     = static_cast<const ITracker&>(tracker);
    //GeomHandle<ITracker> tracker;

    ConditionsHandle<TrackerCalibrations> tcal("ignored");
    const double signalVelocity = tcal->SignalVelocity(StrawIndex(0));//231.;//299.792458; // mm/ns

    // A container to hold the output hits.
    auto_ptr<StrawHitCollection>             strawHits(new StrawHitCollection);
    auto_ptr<StrawHitMCTruthCollection>      truthHits(new StrawHitMCTruthCollection);
    auto_ptr<PtrStepPointMCVectorCollection> mcptrHits(new PtrStepPointMCVectorCollection);

    // Ask the event to give us a handle to the requested hits.
    art::Handle<StepPointMCCollection> points;
    event.getByLabel(_g4ModuleLabel,_trackerStepPoints,points);

    // Organize hits by cells

    //typedef std::map<StepPointMC::VolumeId_type,std::vector<int> > DriftCellHitMap;
    //typedef std::map<StrawIndex,std::vector<int> > DriftCellHitMap;

    DriftCellHitMap hitmap;
    fillHitMap( event, hitmap );


//    for ( size_t i=0; i<points->size(); ++i){
//      StepPointMC const& hit = (*points)[i];
//      if( hit.totalEDep()<_minimumEnergy ) continue; // Skip steps with very low energy deposition
//      //StepPointMC::VolumeId_type dcell_id = hit.volumeId();
//      StrawIndex dcell_id = hit.strawIndex();
//      vector<int> &hits_id = hitmap[dcell_id];
//      hits_id.push_back(i);
//    }

    vector<StepHit> dcell_hits;
    double tdrifterr=0;
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

      //      std::cout<<"cell lay,cell "<<dcell_id.asUint()/10000<<" "<< dcell_id.asUint()%10000<<endl;
      // Prepare info for hit creation
      dcell_hits.clear();

      // Loop over all hits found for this straw

      //vector<int> const& ihits = idcell->second;
      vector<art::Ptr<StepPointMC> > const& ihits = idcell->second;

      for( size_t i=0; i<ihits.size(); i++ ) {

        //int hitRef = ihits[i];
        //StepPointMC const& hit = (*points)[hitRef];
        StepPointMC const& hit = *(ihits.at(i));
        CLHEP::Hep3Vector  const& pos = hit.position();
        CLHEP::Hep3Vector  const& mom = hit.momentum();
        double length  = hit.stepLength();
        double edep    = hit.totalEDep();
        double hitTime = hit.time();

        bool skip_hit;
        //cout<<"hit Z "<<pos.getZ()<<" time "<<hitTime<<endl;
        if (_zHittingZonePresent) {
                int izZone;
                izZone=_nzHittingZone-1;
                if ( pos.getZ()<=_zZoneLimits.at(0) ) {
                        //cout<<"In zone 1"<<endl;
                        if ( _zZoneActive.at(0) ) {
                                if ( hitTime<_zZoneMinLiveTime.at(0) ) { /*cout<<"skipped"<<endl;*/ continue; }
                        } else { /*cout<<"skipped"<<endl;*/ continue; }
                } else if ( pos.getZ()>_zZoneLimits.at(izZone-1) ) {
                        //cout<<"In zone "<<_nzHittingZone<<endl;
                        if ( _zZoneActive.at(izZone) ) {
                                if ( hitTime<_zZoneMinLiveTime.at(izZone) ) { /*cout<<"skipped"<<endl;*/ continue; }
                        } else { /*cout<<"skipped"<<endl;*/ continue; }
                } else {
                        skip_hit = false;
                        for ( izZone=1; izZone<_nzHittingZone-1; izZone++ ) {
                                if ( pos.getZ()>_zZoneLimits.at(izZone-1) && pos.getZ()<=_zZoneLimits.at(izZone) ) {
                                        //cout<<"In zone "<<izZone+1<<endl;
                                        if ( _zZoneActive.at(izZone) ) {
                                                if ( hitTime<_zZoneMinLiveTime.at(izZone) ) { skip_hit = true; break; }
                                        } else { skip_hit = true; break; }
                                }
                        }
                        if (skip_hit) { /*cout<<"skipped"<<endl;*/ continue; }
                }

        }

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
	    //double hit_dca_l = pca.dca();
            hit_pca = pca.point1();
            sign=w.cross(mom).dot(hit_pca-p2);
	    
	    //small correction (~10mum) forr track is not line
	    double dl=(p2-pos).perp();
	    CLHEP::Hep3Vector b(0,0,1.);
	    CLHEP::Hep3Vector pcurv=p2+b.cross(mom)/mom.perp()*dl*dl/2/(mom.perp()*3);
	    hit_dca=(pcurv-hit_pca).mag();
	    // cout<<"hit "<<i<<" dca "<<hit_dca_l<<" "<<hit_dca<<" dl "<<dl<<" "<<dl*dl/2/(mom.perp()*3)<<" pos "<<pos<<" mom "<<mom<<" hit "<<p2<<" w "<<hit_pca<<endl;
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

        D2T d2t;
        tcal->DistanceToTime(dcell_id,hit_dca,mom,d2t);
	tdrifterr=d2t._tdrifterr;
        double driftTime = fabs(d2t._tdrift);
        double distanceToMiddle = (hit_pca-mid).dot(w);
        double hit_t1 = hitTime + driftTime + (strawHalfLength-distanceToMiddle)/signalVelocity;
        //double hit_t2 = -9999.9;//t0 + hitTime + driftTime + (strawHalfLength+distanceToMiddle)/signalVelocity;

        dcell_hits.push_back(StepHit(ihits.at(i),edep,hit_dca*(sign>0?1.:-1.),driftTime,distanceToMiddle,hit_t1));

      } // loop over hits

      // Now that we calculated estimated signal time for all G4Steps, we can analyze
      // the time structure and create StrawHits

      // First we need to sort StepHits according to t1 time
      sort(dcell_hits.begin(), dcell_hits.end() );

      if ( ncalls < _maxFullPrint && _diagLevel > 2 ) {
        for( size_t i=0; i<dcell_hits.size(); i++ ) {
          cout << "MakeDriftCellHit: StepHit # (" << dcell_hits[i]._ptr.id() << " " << dcell_hits[i]._ptr.key() << ")"
               << " DCA=" << dcell_hits[i]._dca
               << " driftT=" << dcell_hits[i]._driftTime
               << " distToMid=" << dcell_hits[i]._distanceToMid
               << " t1=" << dcell_hits[i]._t1
               /*<< " t2=" << dcell_hits[i]._t2*/
               << " edep=" << dcell_hits[i]._edep
               << endl;
        }
      }

      // Now loop over all StepHits and create StrawHits as needed

      if( dcell_hits.size()<1 ) continue; // This should never be needed. Added for safety.

      double digi_time   = dcell_hits[0]._t1+_gaussian.fire(0.,tdrifterr);
      double digi_t2     = 0.0;//dcell_hits[0]._t2;
      double digi_edep   = dcell_hits[0]._edep;
      double digi_driftT = dcell_hits[0]._driftTime;
      double digi_toMid  = dcell_hits[0]._distanceToMid;
      double digi_dca    = dcell_hits[0]._dca;
      double deltadigitime;

      PtrStepPointMCVector mcptr;
      //mcptr.push_back( art::Ptr<StepPointMC>( points, dcell_hits[0]._hit_id ) );
      mcptr.push_back( dcell_hits[0]._ptr );


      for( size_t i=1; i<dcell_hits.size(); i++ ) {
        if( (dcell_hits[i]._t1-dcell_hits[i-1]._t1) > _minimumTimeGap ) {
          // The is bit time gap - save current data as a hit...
          strawHits->push_back(StrawHit(/*StrawIndex(*/dcell_id/*)*/,digi_time,digi_t2/*digi_t2-digi_time*/,digi_edep));
          truthHits->push_back(StrawHitMCTruth(digi_driftT,digi_dca,digi_toMid));
          mcptrHits->push_back(mcptr);
          // ...and create new hit
          mcptr.clear();
          //mcptr.push_back( art::Ptr<StepPointMC>(points, dcell_hits[i]._hit_id));
          mcptr.push_back( dcell_hits[i]._ptr);
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
          //mcptr.push_back( art::Ptr<StepPointMC>(points, dcell_hits[i]._hit_id));
          mcptr.push_back( dcell_hits[i]._ptr );
        }
      }
      //deltadigitime=(digi_t2-digi_time)+_gaussian.fire(0.,_distSigma/_timetodist);
      deltadigitime=0.0;
      strawHits->push_back(StrawHit(/*StrawIndex(*/dcell_id/*)*/,digi_time,deltadigitime/*deltadigitime*/,digi_edep));
      truthHits->push_back(StrawHitMCTruth(digi_driftT,digi_dca,digi_toMid));
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
