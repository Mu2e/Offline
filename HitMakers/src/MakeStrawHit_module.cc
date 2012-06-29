//
// An EDProducer Module that reads StepPointMC objects and turns them into
// StrawHit objects.
//
// $Id: MakeStrawHit_module.cc,v 1.17 2012/06/29 21:30:46 genser Exp $
// $Author: genser $
// $Date: 2012/06/29 21:30:46 $
//
// Original author Rob Kutschke. Updated by Ivan Logashenko.
//                               Updated by Hans Wenzel to include sigma in deltat

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "LTrackerGeom/inc/LTracker.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "Mu2eUtilities/inc/LinePointPCA.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/MassCache.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "SeedService/inc/SeedService.hh"

// art includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Other includes.
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>

using namespace std;

namespace mu2e {

  // Utility class (structure) to hold calculated drift time for G4 hits

  class StepHit {

  public:

    art::Ptr<StepPointMC> _ptr;
    double _edep;
    double _dca;
    double _driftTime;
    double _distanceToMid;
    double _t1;
    double _t2;

    StepHit( art::Ptr<StepPointMC> const& ptr, double edep, double dca, double driftT, double toMid, double t1, double t2):
      _ptr(ptr), _edep(edep), _dca(dca), _driftTime(driftT),
      _distanceToMid(toMid), _t1(t1), _t2(t2) {}

    // This operator is overloaded in order to time-sort the hits
    bool operator <(const StepHit& b) const { return (_t1 < b._t1); }

  };

  //--------------------------------------------------------------------
  //
  //
  class MakeStrawHit : public art::EDProducer {

  public:
    explicit MakeStrawHit(fhicl::ParameterSet const& pset) :

      // Parameters
      _diagLevel(pset.get<int>("diagLevel",0)),
      _maxFullPrint(pset.get<int>("maxFullPrint",5)),
      _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
      _minimumEnergy(pset.get<double>("minimumEnergy",0.0001)), // MeV
      _minimumAmplitude(pset.get<double>("minimumAmplitude",0.0001)), // mV
      _addCrossTalkHits(pset.get<bool>("addCrossTalkHits",false)),
      _minimumLength(pset.get<double>("minimumLength",0.01)),   // mm
      _minimumTimeGap(pset.get<double>("minimumTimeGap",100.0)),// ns
      _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
      _enableFlightTimeCorrection(pset.get<bool>("flightTimeCorrection",false)),

      // Random number distributions
      _gaussian( createEngine( art::ServiceHandle<SeedService>()->getSeed() ) ),

      _messageCategory("HITS"),

      // Control some information messages.
      _firstEvent(true){

      // Tell the framework what we make.
      produces<StrawHitCollection>();
      produces<StrawHitMCTruthCollection>();
      produces<PtrStepPointMCVectorCollection>("StrawHitMCPtr");

    }

    virtual ~MakeStrawHit() { }

    virtual void beginJob();

    void produce( art::Event& e);

  private:

    typedef std::map<StrawIndex,std::vector<art::Ptr<StepPointMC> > > StrawStepPointMap;

    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Name of the tracker StepPoint collection
    std::string _trackerStepPoints;

    // Parameters
    double _minimumEnergy;  // minimum energy deposition of G4 step
    double _minimumAmplitude; // actually min amplitude difference for crosstalk
    bool   _addCrossTalkHits; // should we add cross talk hits?
    double _minimumLength;  // is G4Step is shorter than this, consider it a point
    double _minimumTimeGap;
    string _g4ModuleLabel;  // Name of the module that made these hits.
    bool   _enableFlightTimeCorrection;

    // Random number distributions
    CLHEP::RandGaussQ _gaussian;

    // A category for the error logger.
    const std::string _messageCategory;

    // Give some informationation messages only on the first event.
    bool _firstEvent;

    void fillHitMap ( art::Event const& event, StrawStepPointMap& hitmap );

  };

  void MakeStrawHit::beginJob(){
  }


  // Find StepPointMCs in the event and use them to fill the hit map.
  void MakeStrawHit::fillHitMap ( art::Event const& event, StrawStepPointMap& hitmap){

    // This selector will select only data products with the given instance name.
    art::ProductInstanceNameSelector selector("tracker");

    typedef std::vector< art::Handle<StepPointMCCollection> > HandleVector;

    // Get all of the tracker StepPointMC collections from the event:
    HandleVector stepsHandles;
    event.getMany( selector, stepsHandles);

    // Informational message on the first event.
    if ( _firstEvent ) {
      mf::LogInfo log(_messageCategory);
      log << "MakeStrawHit::fillHitMap will use StepPointMCs from: \n";
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
  MakeStrawHit::produce(art::Event& event) {

    if ( _diagLevel > 1 ) cout << "MakeStrawHit: produce() begin; event " << event.id().event() << endl;

    static int ncalls(0);
    ++ncalls;

    // Get a reference to one of the L or T trackers.
    // Throw exception if not successful.
    const Tracker& tracker = getTrackerOrThrow();

    // Containers to hold the output information.
    auto_ptr<StrawHitCollection>             strawHits(new StrawHitCollection);
    auto_ptr<StrawHitMCTruthCollection>      truthHits(new StrawHitMCTruthCollection);
    auto_ptr<PtrStepPointMCVectorCollection> mcptrHits(new PtrStepPointMCVectorCollection);

    // Handle to the conditions service
    ConditionsHandle<TrackerCalibrations> trackerCalibrations("ignored");

    // Organize the StepPointMCs by straws.
    StrawStepPointMap hitmap;
    fillHitMap( event, hitmap );

    if (_addCrossTalkHits) {

      //
      // loop over all straws and generate crosstalk hits
      //

      StrawStepPointMap hitsToBeAdded;

      std::deque<Straw> const & allStraws = tracker.getAllStraws();

      for ( std::deque<Straw>::const_iterator i = allStraws.begin(), e = allStraws.end(); i!=e; ++i){

        Straw const & straw = *i;

        // loop over all straws; get their nearest nighbours, 
        // collect new hits (StepPointMCs) based on those

        StrawIndex si = straw.index();
        std::vector<StrawIndex> const & nearestStrawInds = straw.nearestNeighboursByIndex();

        if ( ncalls < _maxFullPrint && _diagLevel > 2 ) {
          cout << "MakeStrawHit: Straw #" << si << endl;
        }

        for( size_t nsi=0, nse=nearestStrawInds.size(); nsi!=nse; ++nsi ) {

          if ( ncalls < _maxFullPrint && _diagLevel > 2 ) {
            cout << "MakeStrawHit    neighbour: " << nearestStrawInds[nsi];
          } // if

          StrawStepPointMap::const_iterator neighbourStrawIterator = 
            hitmap.find(nearestStrawInds[nsi]);

          if (neighbourStrawIterator!=hitmap.end()) { 
            // this neighbour straw has hits
            std::vector<art::Ptr<StepPointMC> > const & neighbourStrawHitPtrs = 
              neighbourStrawIterator->second;
            
            if ( ncalls < _maxFullPrint && _diagLevel > 3 ) {
              cout << " has " << neighbourStrawHitPtrs.size() 
                   << " StepPointMC hits; The ones above threshold are: ";
            } // if

            double digi_ampl;
            E2A e2a;

            for( size_t nii=0, eenii=neighbourStrawHitPtrs.size(); nii!=eenii; ++nii ) {

              art::Ptr<StepPointMC> const & thisNeighbHitPtr = neighbourStrawHitPtrs[nii];

              // collect the hit (StepPointMC) if it is above threshold

              // we translate the energy to amplitude here; FIXME this
              // may be slightly "out of step" or oversimplified as we
              // have not formed full hits yet

              trackerCalibrations->EnergyToAmplitude(thisNeighbHitPtr->strawIndex(), 
                                                     thisNeighbHitPtr->totalEDep(), e2a);
              digi_ampl =  (e2a._ampl + _gaussian.fire(0.,e2a._amplerr))*
                trackerCalibrations->CrossTalk(si,thisNeighbHitPtr->strawIndex());

              if (digi_ampl<_minimumAmplitude) continue;

              if ( ncalls < _maxFullPrint && _diagLevel > 3 ) {
                cout << " (" 
                     << neighbourStrawHitPtrs[nii].id() << " " 
                     << neighbourStrawHitPtrs[nii].key() << ")";
              } // if

              hitsToBeAdded[si].push_back(thisNeighbHitPtr);

            } // for
            
          } // if
          if ( ncalls < _maxFullPrint && _diagLevel > 2 ) {
            cout << endl;
          } // if

        } // for

      } // for

      // now actually merge in the crosstalk hits

      // well, but how to calculate t1 etc...? since we have to do it using another straw as its basis...
      // but it may just work if we use proper straw index

      for(StrawStepPointMap::const_iterator istraw = hitsToBeAdded.begin(), istrawe = hitsToBeAdded.end(); 
          istraw != istrawe; ++istraw) {

        // The StepPointMCs to be added that are on this straw.
        std::vector<art::Ptr<StepPointMC> > const& ihits = istraw->second;

        for( std::vector<art::Ptr<StepPointMC> >::const_iterator spp=ihits.begin(), e=ihits.end();
             spp!=e ; ++spp){
          hitmap[istraw->first].push_back(*spp);

        } // for

      } // for

    } // if 

    // Temporary working space per straw.
    std::vector<StepHit> straw_hits;

    // Cache of recently used masses from the particle data table.
    MassCache cache;

    // Loop over all straws and create StrawHits. There can be several
    // hits per straw if they are separated by time. The general algorithm
    // is as follows: calculate signal time for each G4step, order them
    // in time and look for gaps. If gap exceeds _minimumTimeGap create
    // separate hit.

    // note that the crosstalk hits will have different straw index in their StepPointMCs

    // we still need to rescale the hit energy due to the crosstalk; 
    // FIXME we do it at the energy not amplitude level for now

    for(StrawStepPointMap::const_iterator istraw = hitmap.begin(), istrawe = hitmap.end(); 
        istraw != istrawe; ++istraw) {

      if ( ncalls < _maxFullPrint && _diagLevel > 2 ) {
        cout << "MakeStrawHit: straw ID=" << istraw->first 
             << ", " << tracker.getStraw(istraw->first).id()
             << ": number of G4 step hits " << istraw->second.size()
             << endl;
      }

      // Get the straw information, also by reference.
      StrawIndex straw_id = istraw->first;

      // we need to get the straw info from the StepPointMCs as they
      // may be different from the straw being worked on

      // The StepPointMCs that are on this straw.
      std::vector<art::Ptr<StepPointMC> > const& ihits = istraw->second;

      // Prepare working space.
      straw_hits.clear();
      straw_hits.reserve(ihits.size());

      // Loop over all hits found for this straw
      int nn(0);
      for( std::vector<art::Ptr<StepPointMC> >::const_iterator i=ihits.begin(), e=ihits.end();
           i !=e ; ++i, ++nn ){

        StepPointMC const& hit = **i;

        StrawIndex spmcStraw_id      = hit.strawIndex();
        Straw const&  straw          = tracker.getStraw(spmcStraw_id);
        CLHEP::Hep3Vector const& mid = straw.getMidPoint();
        CLHEP::Hep3Vector const& w   = straw.getDirection();
        double strawHalfLength       = straw.getHalfLength();
        double signalVelocity        = trackerCalibrations->SignalVelocity(spmcStraw_id);

        CLHEP::Hep3Vector  const& pos = hit.position();
        CLHEP::Hep3Vector  const& mom = hit.momentum();
        double length  = hit.stepLength();
        double edep    = (straw_id == spmcStraw_id) ?
          hit.totalEDep() :
          hit.totalEDep()*trackerCalibrations->CrossTalk(straw_id,spmcStraw_id);
        double hitTime = hit.time();

        if ( ncalls < _maxFullPrint && _diagLevel > 3 ) {

          cout << "MakeStrawHit: Hit #" << nn << " : length=" << length
               << " energy=" << edep << " time=" << hitTime;
          if (straw_id != spmcStraw_id) {
            cout << " is a crosstalk hit of edep " << hit.totalEDep();
          }
          cout    << endl;
          // we'll print the track info
          cout << "MakeStrawHit: hit.simParticle() id, pdgid, parent id, pdgid : " 
               << hit.simParticle()->id() << ", " << hit.simParticle()->pdgId();
          if ( hit.simParticle()->isSecondary() ) {
            cout << " is a secondary of " << (hit.simParticle()->parent())->id() 
                 << ", " << (hit.simParticle()->parent())->pdgId() << endl;
          } else {
            cout << " is a primary" << endl;
          }
        }

        // Calculate the drift distance and the point on the wire at the end of the dca vector.
        double hit_dca;
        CLHEP::Hep3Vector hit_pca;

        // Length along the step from the start to the dca.
        double hit_s(0.);

        if( length < _minimumLength ) {

          // If step length is very small, consider it a point

          LinePointPCA pca(mid, w, pos);
          hit_dca = pca.dca();
          hit_pca = pca.pca();
          hit_s   = 0.;

        } else {

          // Step is not a point. Calculate the distance between two lines.

          TwoLinePCA pca( mid, w, pos, mom);

          if ( pca.s2() >=0 && pca.s2() <= length ){

            // If the point of closest approach is within the step and wire - thats it.
            hit_dca = pca.dca();
            hit_pca = pca.point1();
            hit_s   = pca.s2();

          } else {


            // The point of closest approach is not within the step. In this case
            // the closes distance should be calculated from the ends

            LinePointPCA pca1(mid, w, pos);
            LinePointPCA pca2(mid, w, pos+mom.unit()*length);
            if( pca1.dca() < pca2.dca() ) {
              hit_dca = pca1.dca();
              hit_pca = pca1.pca();
              hit_s   = 0.;

            } else {
              hit_dca = pca2.dca();
              hit_pca = pca2.pca();
              hit_s   = length;
            }

          }

        } // drift distance calculation

        // Flight time of particle from start of step to the DCA.
        double flightTime = 0.;
        if ( _enableFlightTimeCorrection ){
          double mass = cache.mass( hit.simParticle()->pdgId() );
          double p    = mom.mag();
          double e    = sqrt( p*p + mass*mass);
          double beta = p/e;
          flightTime  = ( beta > 0 ) ? hit_s/beta/CLHEP::c_light : 0.;
        }

        // Calculate signal time. It is Geant4 time + signal propagation time
        // t1 is signal time at positive end (along w vector),
        // t2 - at negative end (opposite to w vector)

	D2T d2t;
	trackerCalibrations->DistanceToTime(spmcStraw_id,hit_dca,mom,d2t);
	// smear the time to account for dispersion and measurement error.  Truncate at 0
        double driftTime = std::max(0.0,d2t._tdrift + _gaussian.fire(0.,d2t._tdrifterr));
        double distanceToMiddle = (hit_pca-mid).dot(w);
        // The convention is that the principle time measurement (t1) corresponds to a measurement
        // at the end of the wire as signed by the wire direction vector. t2 is at the near end.
        double hit_t1 = hitTime + flightTime + driftTime + (strawHalfLength-distanceToMiddle)/signalVelocity;
        double hit_t2 = hitTime + flightTime + driftTime + (strawHalfLength+distanceToMiddle)/signalVelocity;

        straw_hits.push_back( StepHit( *i,edep,hit_dca,driftTime,distanceToMiddle,hit_t1,hit_t2));

      } // loop over hits

      // Now that we calculated estimated signal time for all G4Steps, we can analyze
      // the time structure and create StrawHits

      // First we need to sort StepHits according to t1 time
      sort(straw_hits.begin(), straw_hits.end() );

      if ( ncalls < _maxFullPrint && _diagLevel > 1 ) {
        for( size_t i=0, e=straw_hits.size(); i!=e ;++i ) {
          cout << "MakeStrawHit: StepHit # (" << straw_hits[i]._ptr.id() << " " << straw_hits[i]._ptr.key() << ")"
               << " DCA=" << straw_hits[i]._dca
               << " driftT=" << straw_hits[i]._driftTime
               << " distToMid=" << straw_hits[i]._distanceToMid
               << " t1=" << straw_hits[i]._t1
               << " t2=" << straw_hits[i]._t2
               << " edep=" << straw_hits[i]._edep
               << endl;
        }
      }

      // Now loop over all StepHits and create StrawHits as needed

      if( straw_hits.size()<1 ) continue; // This should never be needed. Added for safety.

      double digi_time   = straw_hits[0]._t1;
      double digi_t2     = straw_hits[0]._t2;
      double digi_edep   = straw_hits[0]._edep;
      double digi_driftT = straw_hits[0]._driftTime;
      // the straw which defined digi_toMid is the same one which defined t1
      double digi_toMid  = straw_hits[0]._distanceToMid;
      double digi_dca    = straw_hits[0]._dca;
      // we prepare ourselves to record which ptr defined t1,t2
      size_t digi_t1i    = 0;
      art::Ptr<StepPointMC> digi_t1ptr = straw_hits[0]._ptr;
      size_t digi_t2i    = 0;
      art::Ptr<StepPointMC> digi_t2ptr = straw_hits[0]._ptr;
      double deltadigitime;
      double distSigma;

      // to be introduced later
      // double digi_ampl;
      // E2A e2a;

      PtrStepPointMCVector mcptr;
      mcptr.push_back( straw_hits[0]._ptr );


      if ( ncalls < _maxFullPrint && _diagLevel > 3 ) {

        cout << "MakeStrawHit: processing StepHit # (" << straw_hits[0]._ptr.id() << 
          " " << straw_hits[0]._ptr.key() << ")"
             << " DCA=" << straw_hits[0]._dca
             << " driftT=" << straw_hits[0]._driftTime
             << " distToMid=" << straw_hits[0]._distanceToMid
             << " t1=" << straw_hits[0]._t1
             << " t2=" << straw_hits[0]._t2
             << " edep=" << straw_hits[0]._edep
             << endl;

        cout << "MakeStrawHit: from straw_hits straw_hits[0]._ptr->simParticle() id, pdgid, parent id, pdgid : " 
             << straw_hits[0]._ptr->simParticle()->id() << ", " << straw_hits[0]._ptr->simParticle()->pdgId();
        if ( straw_hits[0]._ptr->simParticle()->isSecondary() ) {
          cout << " is a secondary of " << (straw_hits[0]._ptr->simParticle()->parent())->id() 
               << ", " << (straw_hits[0]._ptr->simParticle()->parent())->pdgId() << endl;
        } else {
          cout << " is a primary" << endl;
        }

      }

      for( size_t i=1, e=straw_hits.size(); i!=e; ++i ) {
        if( (straw_hits[i]._t1-straw_hits[i-1]._t1) > _minimumTimeGap ) {
          // The is bit time gap - save current data as a hit...

          StrawIndex spmcStraw_id      = straw_hits[i]._ptr->strawIndex();
          Straw const&  straw          = tracker.getStraw(spmcStraw_id);

          double strawHalfLength       = straw.getHalfLength();
          double signalVelocity        = trackerCalibrations->SignalVelocity(spmcStraw_id);

          distSigma = trackerCalibrations->TimeDivisionResolution(spmcStraw_id, (strawHalfLength-digi_toMid)/(2.*strawHalfLength) );
          deltadigitime = (digi_t2-digi_time)+_gaussian.fire(0.,2.*distSigma/signalVelocity);

          // trackerCalibrations->EnergyToAmplitude(spmcStraw_id, digi_edep, e2a);
          // digi_ampl =   e2a._ampl + _gaussian.fire(0.,e2a._amplerr); // this should replace digi_edep

          // we need to record straw_hits[i]._ptr of the "defining" hits
          // the time is defined by the very first stepHit at index 0

          // how about t2 and energy?
          // is ancestry important?

          // how about the xtalk?
          // edep to signal is first though

          // lets start with t2: t2ptr

          if ( ncalls < _maxFullPrint && _diagLevel > 3 ) {

            cout << "MakeStrawHit: StepHit # will go to the next hit(" << straw_hits[i]._ptr.id() << 
              " " << straw_hits[i]._ptr.key() << ")"
                 << " DCA=" << straw_hits[i]._dca
                 << " driftT=" << straw_hits[i]._driftTime
                 << " distToMid=" << straw_hits[i]._distanceToMid
                 << " t1=" << straw_hits[i]._t1
                 << " t2=" << straw_hits[i]._t2
                 << " edep=" << straw_hits[i]._edep
                 << endl;

          }

          // strawHits->push_back(StrawHit(straw_id,digi_time,deltadigitime,digi_ampl));
          strawHits->push_back(StrawHit(straw_id,digi_time,deltadigitime,digi_edep));
          truthHits->push_back(StrawHitMCTruth(digi_driftT,digi_dca,digi_toMid));
          mcptrHits->push_back(mcptr);

          if ( ncalls < _maxFullPrint && _diagLevel > 2 ) {

            // will print just added hit
            cout << "MakeStrawHit: added hit: Straw #" << ((*strawHits).back()).strawIndex()
                 << " time="  << ((*strawHits).back()).time()
                 << " dt="    << ((*strawHits).back()).dt()
                 << " edep="  << ((*strawHits).back()).energyDep()
                 << " nsteps="<< ((*mcptrHits).back()).size()
                 << " t1i="   << digi_t1i
              //   << " t1ptr= "<< digi_t1ptr
                 << " t2i="   << digi_t2i
              //   << " t2ptr= "<< digi_t2ptr
                 << endl;

          }

          // ...and create new hit
          mcptr.clear();
          mcptr.push_back( straw_hits[i]._ptr );
          digi_time   = straw_hits[i]._t1;
          digi_t2     = straw_hits[i]._t2;
          digi_edep   = straw_hits[i]._edep;
          digi_driftT = straw_hits[i]._driftTime;
          digi_toMid  = straw_hits[i]._distanceToMid;
          digi_dca    = straw_hits[i]._dca;
          digi_t1i    = i;
          digi_t1ptr  = straw_hits[i]._ptr;
          digi_t2i    = i;
          digi_t2ptr  = straw_hits[i]._ptr;
        } else {
          // Append existing hit

          if ( ncalls < _maxFullPrint && _diagLevel > 3 ) {

            cout << "MakeStrawHit: appending StepHit # (" << straw_hits[i]._ptr.id() << 
              " " << straw_hits[i]._ptr.key() << ")"
                 << " DCA=" << straw_hits[i]._dca
                 << " driftT=" << straw_hits[i]._driftTime
                 << " distToMid=" << straw_hits[i]._distanceToMid
                 << " t1=" << straw_hits[i]._t1
                 << " t2=" << straw_hits[i]._t2
                 << " edep=" << straw_hits[i]._edep
                 << endl;

            cout << "MakeStrawHit: from straw_hits straw_hits[i]._ptr->simParticle() id, pdgid, parent id, pdgid : " 
                 << straw_hits[i]._ptr->simParticle()->id() << ", " << straw_hits[i]._ptr->simParticle()->pdgId();
            if ( straw_hits[i]._ptr->simParticle()->isSecondary() ) {
              cout << " is a secondary of " << (straw_hits[i]._ptr->simParticle()->parent())->id() 
                   << ", " << (straw_hits[i]._ptr->simParticle()->parent())->pdgId() << endl;
            } else {
              cout << " is a primary" << endl;
            }

          }

          if( digi_t2 > straw_hits[i]._t2 ) {
            digi_t2    = straw_hits[i]._t2;
            digi_t2i   = i;
            digi_t2ptr = straw_hits[i]._ptr;
          }

          digi_edep += straw_hits[i]._edep;
          mcptr.push_back( straw_hits[i]._ptr );
        }
      }

      // last hit

      StrawIndex spmcStraw_id      = mcptr.back()->strawIndex();
      Straw const&  straw          = tracker.getStraw(spmcStraw_id);

      double strawHalfLength       = straw.getHalfLength();
      double signalVelocity        = trackerCalibrations->SignalVelocity(spmcStraw_id);

      distSigma = trackerCalibrations->TimeDivisionResolution(spmcStraw_id, (strawHalfLength-digi_toMid)/(2.*strawHalfLength) );
      deltadigitime = (digi_t2 - digi_time)+_gaussian.fire(0.,2.*distSigma/signalVelocity);
      
      // trackerCalibrations->EnergyToAmplitude(spmcStraw_id, digi_edep, e2a);
      // digi_ampl =   e2a._ampl + _gaussian.fire(0.,e2a._amplerr); // this should replace digi_edep

      // strawHits->push_back(StrawHit(straw_id,digi_time,deltadigitime,digi_ampl));
      strawHits->push_back(StrawHit(straw_id,digi_time,deltadigitime,digi_edep));
      // truthHits->push_back(StrawHitMCTruth(digi_driftT,digi_dca,digi_toMid,
      //                      digi_edep,digi_t1ptr,digi_t2ptr));
      truthHits->push_back(StrawHitMCTruth(digi_driftT,digi_dca,digi_toMid));
      mcptrHits->push_back(mcptr);

      if ( ncalls < _maxFullPrint && _diagLevel > 2 ) {

        // will print just added hit
        cout << "MakeStrawHit: added hit: Straw #" << ((*strawHits).back()).strawIndex()
             << " time="  << ((*strawHits).back()).time()
             << " dt="    << ((*strawHits).back()).dt()
             << " edep="  << ((*strawHits).back()).energyDep()
             << " nsteps="<< ((*mcptrHits).back()).size()
             << " t1i="   << digi_t1i
          //  << " t1ptr= "<< digi_t1ptr
             << " t2i="   << digi_t2i
          //  << " t2ptr= "<< digi_t2ptr
             << endl;
        
      }

    }

    if ( ncalls < _maxFullPrint && _diagLevel > 1 ) {
      cout << "MakeStrawHit: Total number of straw hits = " << strawHits->size() << endl;

      for( size_t i=0, e=strawHits->size(); i!=e; ++i ) {
        cout << "MakeStrawHit: Straw #" << (*strawHits)[i].strawIndex()
             << " time="  << (*strawHits)[i].time()
             << " dt="    << (*strawHits)[i].dt()
             << " edep="  << (*strawHits)[i].energyDep()
             << " nsteps="<< (*mcptrHits)[i].size()
             << endl;
      }
    }

    // we have made many changes; we need to check if the sizes are the same

    assert (strawHits->size()==truthHits->size());
    assert (strawHits->size()==mcptrHits->size());

    event.put(strawHits);
    event.put(truthHits);
    event.put(mcptrHits,"StrawHitMCPtr");

    if ( _diagLevel > 1 ) cout << "MakeStrawHit: produce() end" << endl;

    // Done with the first event; disable some messages.
    _firstEvent = false;

  } // end MakeStrawHit::produce.

} // end namespace mu2e

using mu2e::MakeStrawHit;
DEFINE_ART_MODULE(MakeStrawHit)
