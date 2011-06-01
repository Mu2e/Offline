#ifndef EarlyPatRecPR_Incident_hh
#define EarlyPatRecPR_Incident_hh

//=============================================================================
//
// An Incident is a collection of StawHitSynopsis corresponding to some
// limited time window within the event.  
//
// $Id: Incident.hh,v 1.1 2011/06/01 23:29:47 mf Exp $
// $Author: mf $
// $Date: 2011/06/01 23:29:47 $
//
// Original author: Mark Fischler
//
//=============================================================================

// C++ includes
#include <iostream>
#include <vector>

// Framework includes.

// Mu2e includes.
// #include "TrackerGeom/inc/StrawHit.hh"
#include "TrackerGeom/inc/Tracker.hh"

#include "EarlyPatRec/inc/StrawHitSynopsis.hh"

namespace mu2e {

  // Incident:  
  //
  // An Incident is a collection of StrawHits, all in one event
  // with timings such that their time of the pulse reaching the wire
  // lies, for all included hits, in a given time window.  The idea is
  // that this window is chosen such that all the hits are consistent 
  // with being part of the same track.  

  class Incident 
  {
  public:
    // ctor -- provide the TTracker geometry, the hits, and the time range
    Incident(Tracker const & tracker, 
             std::vector<StrawHit const *> & strawHitPtr,
             double t1, double t2 );
    
    Incident() : hits(), _tt(0) {}

    // main use of this class:  
    std:: vector<StrawHitSynopsis> const & getHits() const;
        
  private:
    // work methods
    bool   strawHitIsInTimeWindow( StrawHit const& h,
                                   double t1, double t2 ) const;  
    double rough_t0_star ( StrawHit const& str ) const;
    void   populateDetectorEventModel();
       
    // private data
    std:: vector<StrawHitSynopsis> hits;
    Tracker const * _tt;
    
  }; // Incident class definition

} // end of mu2e namespace

#endif // EarlyPatRecPR_Incident_hh
