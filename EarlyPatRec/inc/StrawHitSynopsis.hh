#ifndef EarlyPatRecPR_StrawHitSynopsis_hh
#define EarlyPatRecPR_StrawHitSynopsis_hh

//=============================================================================
//
// A StawHitSynopsis provides qick information about the hit, plus a pointer
// to the actual straw hit.
//
// $Id: StrawHitSynopsis.hh,v 1.1 2011/06/01 23:29:47 mf Exp $
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
#include "TrackerGeom/inc/Tracker.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"


namespace mu2e {

  class StrawHitSynopsis 
  {
  public:
    // ctor -- provide the TTracker geometry, the hit, and time estimate
    StrawHitSynopsis(Tracker const & tracker, 
             StrawHit const& h,
             double t0_star)
      : _tt(&tracker)
      , _strawHitPtr (&h)
      , _t0_star (t0_star) {}
      
    // TODO - There will be more content to this class; also accessors
    //        (which is why _tt will be needed) 

   StrawHitSynopsis() 
   : _tt(0)
   , _strawHitPtr (0) 
   , _t0_star (0) {}
   
  private:
   // private data
    Tracker const * _tt;
    StrawHit const * _strawHitPtr;
    double _t0_star;
    
  }; // StrawHitSynopsis class definition

} // end of namespace mu2e

#endif // EarlyPatRecPR_StrawHitSynopsis_hh
