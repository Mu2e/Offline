//=============================================================================
//
// Code to produce a vector<SubEvent> from the collection of StrawHits
// for an event.
//
// $Id: SubEventsMaker_module.cc,v 1.6 2013/03/15 15:52:03 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:03 $
//
// Original author: Mark Fischler
//
//=============================================================================

// C++ includes
#include <vector>
#include <cassert>

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


// Root includes.

// Mu2e includes.
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "Mu2eUtilities/inc/SimParticlesWithHits.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "TrackerGeom/inc/StrawId.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "RecoDataProducts/inc/SubEvent.hh"
#include "RecoDataProducts/inc/SubEventCollection.hh"

// TODO -- temporary printout define; replace with actual controls
// #define DIAGNOSICS_FOR_SEM_TERSE



namespace mu2e {;

  // SubEventsMaker -- a producer putting a vector of SubEvents into
  //                   the event.
  //
  // An SubEvent is a collection of StrawHitSynopsis, all in one event
  // with timings such that they are all consistent with being part
  // of the same track.  StrawHitSynopsis is fundamentally a StrawHit
  // (kept by pointer) with some additional information of convenience
  // for early pattern recognition purposes.
  //
  // The vector of SubEvents for an event can be analyzed
  // one entry at a time, to get all the plausible tracks.
  //
  // Depending on how we separate the hits into possible time windows,
  // a hit may appear in more than one SubEvent in the vector, and 
  // it is conceivable that the same track will appear multiple times
  // in multiple SubEvents.  This will need to be sorted out at a later
  // stage in reconstruction, by elimination of near-duplicate tracks.

  class SubEventsMaker : public art::EDProducer
  {
  public:
    explicit SubEventsMaker(fhicl::ParameterSet const& pSet);
    virtual ~SubEventsMaker() { }
    virtual void produce(art::Event& e );
    
  private:

    // helper functions
    void fillSubEvents( art::Event& event,
                        TTracker const& tt,
                        SubEventCollection & SubEvents );
    void fillControlFromPset(fhicl::ParameterSet const& pSet);
    
    // helper sub-class
    class SubEventSeparator 
    {
    public: 
      SubEventSeparator (TTracker const& trk, StrawHitCollection const & hits);
      void timeSlice (double t1, double t2, SubEvent & subEvent) const;
      double rough_t0_star(StrawHit const& h) const;
    private:
      TTracker const& tt;
      StrawHitCollection const & h;
      std::vector<double> t0_star;
    };

    // private member data
    double _timeWindow;
    double _timeStart;
    double _timeEnd;
    double _timeBinStep;
    unsigned int _minSubEventHits;
    std::string  _makerModuleLabel;
   
  }; // SubEventsMaker class definition
 
  // Methods that happen at construction time
  // -----------------------------------------
  
  SubEventsMaker::SubEventsMaker(fhicl::ParameterSet const& pSet)
  {
    produces<SubEventCollection>();
    fillControlFromPset (pSet);     
  }

  void SubEventsMaker::fillControlFromPset(fhicl::ParameterSet const& pSet)
  {
    _timeWindow    = pSet.get<double>("timeWindow",  150.0);
    _timeStart     = pSet.get<double>("timeStart",   675.0);
    _timeEnd       = pSet.get<double>("timeEnd",    1695.0);
    _timeBinStep   = pSet.get<double>("timeBinStep",  50.0);
    _minSubEventHits  =  pSet.get<unsigned int> ("minSubEventHits",6); 
    _makerModuleLabel =  pSet.get<std::string>  ("makerModuleLabel"); 
    
  } // SubEventsMaker::fillControlFromPset(pSet)

  // Methods that happen for each event
  // ----------------------------------
  
  void SubEventsMaker::produce(art::Event& event) 
  {
    GeomHandle<TTracker> ttracker;
    TTracker const& trk(*ttracker);
    std::unique_ptr <SubEventCollection> p (new SubEventCollection);
    SubEventCollection & subEvents(*p);
    fillSubEvents(event, trk, subEvents);
    event.put(std::move(p));
  }
    
  void SubEventsMaker::fillSubEvents( art::Event& event,
                                      TTracker const& trk,
                                      SubEventCollection & subEvents )
  {
    art::Handle<StrawHitCollection> pdataHandle;
    event.getByLabel(_makerModuleLabel, pdataHandle);
    StrawHitCollection const* hitsPtr = pdataHandle.product();
    StrawHitCollection const & hits(*hitsPtr);

#ifdef DIAGNOSICS_FOR_SEM_TERSE
    std::cout << "Filling SubEvents for an event with " 
              << hits.size() << " hits\n";
#endif

    SubEventSeparator separator ( trk, hits );

    for ( double t1 = _timeStart; 
          t1 + _timeWindow - _timeBinStep <= _timeEnd;
          t1 +=  _timeBinStep  )
      {
        double t2 = std::min (t1 + _timeWindow, _timeEnd);
        SubEvent subEvent(t1,t2);
        separator.timeSlice(t1,t2,subEvent);
        if ( subEvent.nHits() >= _minSubEventHits )
          {
            // TODO - this is a temporary diagnostic output
            // We should perhaps add a diagLevel control over this.
#ifdef DIAGNOSICS_FOR_SEM
            std::cout << "About to push back a SubEvent with "
                      << subEvent.nHits() << " hits\n";
#endif 
            subEvents.push_back ( subEvent ); 
#ifdef DIAGNOSICS_FOR_SEM
            std::cout << "pushed back a SubEvent with "
                      << subEvent.nHits() << " hits\n";
#endif 
          } 
        // TODO -- this is a temproary diagnostic printout
#ifdef DIAGNOSICS_FOR_SEM_TERSE
        std::cout << "fillSubEvents: SubEvent in time range (" << t1
                  << ", " << t2 << ") contains " << subEvent.nHits() << " hits\n";
#endif 
      } 
    // TODO TWEAK:  Investigate best choices for time windows, or possibly
    //              different algorithms than this simple one for separating 
    //              into windows 

  } // SubEventsMaker::fillSubEvents(...)


  SubEventsMaker::SubEventSeparator::SubEventSeparator 
  (TTracker const& trk, StrawHitCollection const & hits)
    : tt(trk)
    , h(hits)
    , t0_star()
  {
    typedef StrawHitCollection::const_iterator SHCci;
    for (SHCci hit = h.begin(); hit != h.end(); ++hit)
      {
        t0_star.push_back(rough_t0_star(*hit));
      }
  } // ctor of SubEventSeparator

  double
  SubEventsMaker::SubEventSeparator::rough_t0_star(StrawHit const& h) const
  {
    const double c = 299.8; // mm/ns
    const double wirePropogationSpeed = 200.0; // mm/ns
    const double betaC = c;
    const double averageSecTheta = 1.635; // based on avg tan theta of 1.294
    double strawHalfLength = tt.getStraw(h.strawIndex()).getHalfLength();
    double strawZ =  tt.getStraw(h.strawIndex()).getMidPoint().z();
    double t = h.time();
    double delta_t = h.dt();
    // TODO - this is a temporary diagnostic output
#ifdef DIAGNOSICS_FOR_SEM
    std::cout << "hit at time " << t
              << " with dt = " << delta_t
              << " with hL = " << strawHalfLength
              << " and Z = "   << strawZ
              << " has t0_star = " << 
      t - 0.5*delta_t + strawHalfLength/wirePropogationSpeed
      + strawZ*averageSecTheta/betaC << "\n";
#endif 
    return t - 0.5*delta_t + strawHalfLength/wirePropogationSpeed
      + strawZ*averageSecTheta/betaC;
  } // SubEventSeparator::rough_t0_star

  void
  SubEventsMaker::SubEventSeparator::timeSlice 
  ( double t1, double t2, SubEvent & subEvent) const
  {
    assert(h.size() == t0_star.size()); // elements of t0_star correspond
                                        // to the hits elements of h
    typedef StrawHitCollection::const_iterator SHCci;
    std::vector<double>::const_iterator timeList_it(t0_star.begin());
    for (SHCci h_it = h.begin();
         h_it != h.end(); 
         ++h_it, ++timeList_it) 
      {
        double t = *timeList_it;
        if ( (t < t1) || (t > t2) ) continue;
        subEvent.push_back (StrawHitSynopsis(*h_it, t));
      }
  }  // SubEventSeparator::timeSlice 
   
} // end of namespace mu2e

using mu2e::SubEventsMaker;
DEFINE_ART_MODULE(SubEventsMaker);
