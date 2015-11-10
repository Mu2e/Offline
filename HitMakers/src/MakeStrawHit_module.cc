//
// An EDProducer Module that reads StepPointMC objects and turns them into
// StrawHit objects.
//
// $Id: MakeStrawHit_module.cc,v 1.24 2014/05/30 19:15:32 brownd Exp $
// $Author: brownd $
// $Date: 2014/05/30 19:15:32 $
//
// Original author Rob Kutschke. Updated by Ivan Logashenko.
//                               Updated by Hans Wenzel to include sigma in deltat
//                               Updated by K.Genser added x-talk, conversion to amplitude, restructured

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GlobalConstantsService/inc/MassCache.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "SeedService/inc/SeedService.hh"

#include "MCDataProducts/inc/StepPointMCStrawHit.hh"
#include "HitMakers/inc/formStepPointMCStrawHit.hh"
#include "TrackerConditions/inc/DeadStrawList.hh"
#include "GeneralUtilities/inc/TwoLinePCA.hh"

// art includes.
#include "art/Persistency/Common/Ptr.h"
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

  class MakeStrawHit : public art::EDProducer {

  public:
    explicit MakeStrawHit(fhicl::ParameterSet const& pset);
    // Accept compiler written d'tor.

    virtual void beginJob();
    virtual void beginRun( art::Run& run );

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

    // List of dead straws as a parameter set; needed at beginRun time.
    fhicl::ParameterSet _deadStraws;

    // Access the live/dead status of the straw.
    DeadStrawList _strawStatus;

    void fillHitMap ( art::Event const& event, StrawStepPointMap& hitmap );

  };

  MakeStrawHit::MakeStrawHit(fhicl::ParameterSet const& pset) :

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
    _firstEvent(true),

    _deadStraws(pset.get<fhicl::ParameterSet>("deadStrawList", fhicl::ParameterSet())),
    _strawStatus(_deadStraws) {

    // Tell the framework what we make.
    produces<StrawHitCollection>();
    produces<StrawHitMCTruthCollection>();
    produces<PtrStepPointMCVectorCollection>("StrawHitMCPtr");

  }

  void MakeStrawHit::beginJob(){
  }

  void MakeStrawHit::beginRun( art::Run& run ){
    _strawStatus.reset(_deadStraws);
  }


  // Find StepPointMCs in the event and use them to fill the hit map.
  void MakeStrawHit::fillHitMap ( art::Event const& event, StrawStepPointMap& hitmap){

    const Tracker& tracker = getTrackerOrThrow();

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

        // Skip dead straws.
        if ( _strawStatus.isDead(step.strawIndex()) ) continue;

        // Skip steps that occur in the deadened region near the end of each wire.
        Straw const& straw = tracker.getStraw(step.strawIndex());
        TwoLinePCA pca( straw.getMidPoint(), straw.getDirection(), step.position(), step.momentum().unit() );
        double zcut = pca.s1()/straw.getDetail().activeHalfLength();
        if ( std::abs(zcut)>1.0 ) continue;


        // The main event for this method.
        StrawIndex const & straw_id = step.strawIndex();
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
    unique_ptr<StrawHitCollection>             strawHits(new StrawHitCollection);
    unique_ptr<StrawHitMCTruthCollection>      truthHits(new StrawHitMCTruthCollection);
    unique_ptr<PtrStepPointMCVectorCollection> mcptrHits(new PtrStepPointMCVectorCollection);

    // Handle to the conditions service
    ConditionsHandle<TrackerCalibrations> trackerCalibrations("ignored");

    // Organize the StepPointMCs by straws.
    StrawStepPointMap hitmap;
    fillHitMap( event, hitmap );

    if (_addCrossTalkHits) {

      //
      // loop over all straws and add crosstalk StepPointMCs
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

              // we may also need to store the crosstalk value in a "truth" object

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

      // we have to do calculations using the right straw index as its basis...


      for(StrawStepPointMap::const_iterator 
            istraw = hitsToBeAdded.begin(), istrawe = hitsToBeAdded.end(); 
          istraw != istrawe; ++istraw) {

        // The StepPointMCs to be added that are on this straw.
        std::vector<art::Ptr<StepPointMC> > const& ihits = istraw->second;

        for( std::vector<art::Ptr<StepPointMC> >::const_iterator spp=ihits.begin(), e=ihits.end();
             spp!=e ; ++spp){
          hitmap[istraw->first].push_back(*spp);

        } // for

      } // for

    } // if (_addCrossTalkHits)

    // Temporary working space per straw.
    std::vector<StepPointMCStrawHit> spmcStraw_Hits;

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
      StrawIndex const & straw_id = istraw->first;

      // we sometimes need to get the straw info from the StepPointMCs as they
      // may be different from the straw being worked on

      // The StepPointMCs that are on this straw.
      std::vector<art::Ptr<StepPointMC> > const& ihits = istraw->second;

      // Prepare working space.
      spmcStraw_Hits.clear();
      spmcStraw_Hits.reserve(ihits.size());

      // Loop over all StepPointMCs found for this straw
      int nn(0);
      for( std::vector<art::Ptr<StepPointMC> >::const_iterator spmci=ihits.begin(), e=ihits.end();
           spmci !=e ; ++spmci, ++nn ){

        art::Ptr<StepPointMC> const spmcptr = *spmci;
        std::unique_ptr<StepPointMCStrawHit> const tmpshp = 
          formStepPointMCStrawHit(
                                  spmcptr,
                                  straw_id,
                                  _minimumLength,
                                  _enableFlightTimeCorrection,
                                  cache,
                                  _gaussian,
                                  tracker,
                                  trackerCalibrations);
        
        spmcStraw_Hits.push_back(*tmpshp);

        if ( ncalls < _maxFullPrint && _diagLevel > 3 ) {

          StepPointMC const & hit = **spmci;
          StrawIndex const & spmcStraw_id = hit.strawIndex();

          double edep    = (straw_id == spmcStraw_id) ?
            hit.totalEDep() :
            hit.totalEDep()*trackerCalibrations->CrossTalk(straw_id,spmcStraw_id);

          cout << "MakeStrawHit: Hit #" << nn << " : length=" << hit.stepLength()
               << " energy=" << edep << " time=" << hit.time();
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

      } // loop over StepPointMCs

      // Now that we calculated estimated signal time for all G4Steps, we can analyze
      // the time structure and create StrawHits

      // First we need to sort StepPointMCStrawHits according to t1 time
      sort(spmcStraw_Hits.begin(), spmcStraw_Hits.end() );

      if ( ncalls < _maxFullPrint && _diagLevel > 1 ) {
        for( size_t i=0, e=spmcStraw_Hits.size(); i!=e ;++i ) {
          StepPointMCStrawHit const & spmcStrawhit = spmcStraw_Hits[i];
          cout << "MakeStrawHit: StepPointMCStrawHit # (" 
               << spmcStrawhit._ptr.id() << " " << spmcStrawhit._ptr.key() << ")"
               << " DCA=" << spmcStrawhit._dca
               << " driftTNonSm=" << spmcStrawhit._driftTimeNonSm
               << " driftT=" << spmcStrawhit._driftTime
               << " distToMid=" << spmcStrawhit._distanceToMid
               << " t1=" << spmcStrawhit._t1
               << " t2=" << spmcStrawhit._t2
               << " edep=" << spmcStrawhit._edep
               << endl;
        }
      }

      // Now loop over all StepPointMCStrawHits and create StrawHits as needed

      double digi_ampl;
      E2A e2a;

      if( spmcStraw_Hits.size()<1 ) continue; // This should never be needed. Added for safety.

      StepPointMCStrawHit const & spmcStrawHit = spmcStraw_Hits[0];

      double digi_time   = spmcStrawHit._t1;
      double digi_t2     = spmcStrawHit._t2;
      double digi_edep   = spmcStrawHit._edep;
      double digi_driftTNonSm = spmcStrawHit._driftTimeNonSm;
      // will use it some day double digi_driftT = spmcStrawHit._driftTime;
      // the straw which defined digi_toMid is the same one which defined t1
      double digi_toMid  = spmcStrawHit._distanceToMid;
      double digi_dca    = spmcStrawHit._dca;
      // we prepare ourselves to record which ptr defined t1,t2
      size_t digi_t1i    = 0;
      art::Ptr<StepPointMC> digi_t1ptr = spmcStrawHit._ptr;
      size_t digi_t2i    = 0;
      art::Ptr<StepPointMC> digi_t2ptr = spmcStrawHit._ptr;

      double deltadigitime;
      double distSigma;
      double strawHalfLength;
      double signalVelocity;

      PtrStepPointMCVector mcptr;
      mcptr.push_back( spmcStrawHit._ptr );

      if ( ncalls < _maxFullPrint && _diagLevel > 3 ) {

        cout << "MakeStrawHit: processing StepPointMCStrawHit # (" << spmcStrawHit._ptr.id() 
             << " " << spmcStrawHit._ptr.key() << ")"
             << " DCA=" << spmcStrawHit._dca
             << " driftTNonSm=" << spmcStrawHit._driftTimeNonSm
             << " driftT=" << spmcStrawHit._driftTime
             << " distToMid=" << spmcStrawHit._distanceToMid
             << " t1=" << spmcStrawHit._t1
             << " t2=" << spmcStrawHit._t2
             << " edep=" << spmcStrawHit._edep;

        if (straw_id != spmcStrawHit._ptr->strawIndex()) {
          cout << " a crosstalk hit ";
        }
        cout    << endl;

        cout << "MakeStrawHit: from spmcStraw_Hits spmcStraw_Hits[0]._ptr->simParticle() "
             << "id, pdgid, parent id, pdgid : " 
             << spmcStrawHit._ptr->simParticle()->id() 
             << ", " << spmcStrawHit._ptr->simParticle()->pdgId();
        if ( spmcStrawHit._ptr->simParticle()->isSecondary() ) {
          cout << " is a secondary of " << (spmcStrawHit._ptr->simParticle()->parent())->id() 
               << ", " << (spmcStrawHit._ptr->simParticle()->parent())->pdgId() << endl;
        } else {
          cout << " is a primary" << endl;
        }

      }

      for( size_t spmci=1, e=spmcStraw_Hits.size(); spmci!=e; ++spmci ) {

        StepPointMCStrawHit const & spmcStrawHit = spmcStraw_Hits[spmci];
        StrawIndex const & spmcStraw_id = spmcStrawHit._ptr->strawIndex();

        if( (spmcStrawHit._t1-spmcStraw_Hits[spmci-1]._t1) > _minimumTimeGap ) {
          // The is bit time gap - save current data as a hit...

          Straw const & straw = tracker.getStraw(straw_id);

          // FIXME we may need to do average straw lenght, distSigma,
          // signal velocity etc... if t1/t2 was calculated based on
          // two different straws as in case of the crosstalk hits

          strawHalfLength = straw.getHalfLength();
          signalVelocity  = trackerCalibrations->SignalVelocity(straw_id);

          distSigma = trackerCalibrations->
            TimeDivisionResolution(straw_id, (strawHalfLength-digi_toMid)/(2.*strawHalfLength) );
          deltadigitime = (digi_t2-digi_time)+_gaussian.fire(0.,2.*distSigma/signalVelocity);

          // we translate the energy to amplitude after all energies were added
          // we do it for the specific straw (not hit)
          trackerCalibrations->EnergyToAmplitude(straw_id, digi_edep, e2a);
          // this should replace digi_edep
          digi_ampl = e2a._ampl + _gaussian.fire(0.,e2a._amplerr); 

          // note that crosstalk hits energy was adjusted in formStepPointMCStrawHit

          // we "record" spmcStrawHit._ptr of the "defining" hits
          // the time is defined by the very first stepHit at index 0
          // t2 is more complicated

          if ( ncalls < _maxFullPrint && _diagLevel > 3 ) {

            cout << "MakeStrawHit: StepPointMCStrawHit # will go to the next hit(" 
                 << spmcStrawHit._ptr.id() << " "
                 << spmcStrawHit._ptr.key() << ")"
                 << " DCA=" << spmcStrawHit._dca
                 << " driftTNonSm=" << spmcStrawHit._driftTimeNonSm
                 << " driftT=" << spmcStrawHit._driftTime
                 << " distToMid=" << spmcStrawHit._distanceToMid
                 << " t1=" << spmcStrawHit._t1
                 << " t2=" << spmcStrawHit._t2
                 << " edep=" << spmcStrawHit._edep;
            if (straw_id != spmcStraw_id) {
              cout << " a crosstalk hit ";
            }
            cout    << endl;

          }

          strawHits->push_back(StrawHit(straw_id,digi_time,deltadigitime,digi_ampl));
          // strawHits->push_back(StrawHit(straw_id,digi_time,deltadigitime,digi_edep));
          truthHits->push_back(StrawHitMCTruth(digi_driftTNonSm,digi_dca,digi_toMid));
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
          mcptr.push_back( spmcStrawHit._ptr );
          digi_time   = spmcStrawHit._t1;
          digi_t2     = spmcStrawHit._t2;
          digi_edep   = spmcStrawHit._edep;
          digi_driftTNonSm = spmcStrawHit._driftTimeNonSm;
          // will use it some day digi_driftT = spmcStrawHit._driftTime;
          digi_toMid  = spmcStrawHit._distanceToMid;
          digi_dca    = spmcStrawHit._dca;
          digi_t1i    = spmci;
          digi_t1ptr  = spmcStrawHit._ptr;
          digi_t2i    = spmci;
          digi_t2ptr  = spmcStrawHit._ptr;
        } else {
          // Append existing hit

          if( digi_t2 > spmcStrawHit._t2 ) {
            digi_t2    = spmcStrawHit._t2;
            digi_t2i   = spmci;
            digi_t2ptr = spmcStrawHit._ptr;
            // should one not recalculate delatadigitime???

          }

          digi_edep += spmcStrawHit._edep;
          mcptr.push_back( spmcStrawHit._ptr );

          if ( ncalls < _maxFullPrint && _diagLevel > 3 ) {

            cout << "MakeStrawHit: appending StepPointMCStrawHit # (" 
                 << spmcStrawHit._ptr.id() << " "
                 << spmcStrawHit._ptr.key() << ")"
                 << " DCA=" << spmcStrawHit._dca
                 << " driftTNonSm=" << spmcStrawHit._driftTimeNonSm
                 << " driftT=" << spmcStrawHit._driftTime
                 << " distToMid=" << spmcStrawHit._distanceToMid
                 << " t1=" << spmcStrawHit._t1
                 << " t2=" << spmcStrawHit._t2
                 << " edep=" << spmcStrawHit._edep;
            if (straw_id != spmcStraw_id) {
              cout << " a crosstalk hit ";
            }
            cout    << endl;

            cout << "MakeStrawHit: from spmcStraw_Hits spmcStraw_Hits[spmci]._ptr->simParticle() "
                 << "id, pdgid, parent id, pdgid : " 
                 << spmcStrawHit._ptr->simParticle()->id() 
                 << ", " << spmcStrawHit._ptr->simParticle()->pdgId();
            if ( spmcStrawHit._ptr->simParticle()->isSecondary() ) {
              cout << " is a secondary of " << (spmcStrawHit._ptr->simParticle()->parent())->id() 
                   << ", " << (spmcStrawHit._ptr->simParticle()->parent())->pdgId() << endl;
            } else {
              cout << " is a primary" << endl;
            }
          }
        }
      }

      // last hit (calculations here are specific to the given straw)

      Straw const & straw             = tracker.getStraw(straw_id);

      strawHalfLength = straw.getHalfLength();
      signalVelocity  = trackerCalibrations->SignalVelocity(straw_id);

      distSigma = trackerCalibrations->
        TimeDivisionResolution(straw_id, (strawHalfLength-digi_toMid)/(2.*strawHalfLength) );
      deltadigitime = (digi_t2 - digi_time)+_gaussian.fire(0.,2.*distSigma/signalVelocity);
      
      trackerCalibrations->EnergyToAmplitude(straw_id, digi_edep, e2a);
      digi_ampl =   e2a._ampl + _gaussian.fire(0.,e2a._amplerr); 

      strawHits->push_back(StrawHit(straw_id,digi_time,deltadigitime,digi_ampl));
      // strawHits->push_back(StrawHit(straw_id,digi_time,deltadigitime,digi_edep));
      // truthHits->push_back(StrawHitMCTruth(digi_driftTNonSm,digi_dca,digi_toMid,
      //                      digi_edep,digi_t1ptr,digi_t2ptr));
      truthHits->push_back(StrawHitMCTruth(digi_driftTNonSm,digi_dca,digi_toMid));
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

    event.put(std::move(strawHits));
    event.put(std::move(truthHits));
    event.put(std::move(mcptrHits),"StrawHitMCPtr");

    if ( _diagLevel > 1 ) cout << "MakeStrawHit: produce() end" << endl;

    // Done with the first event; disable some messages.
    _firstEvent = false;

  } // end MakeStrawHit::produce.

} // end namespace mu2e

using mu2e::MakeStrawHit;
DEFINE_ART_MODULE(MakeStrawHit)
