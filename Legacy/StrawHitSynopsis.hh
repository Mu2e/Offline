#ifndef EarlyPatRecPR_StrawHitSynopsis_hh
#define EarlyPatRecPR_StrawHitSynopsis_hh

//=============================================================================
//
// A StawHitSynopsis provides qick information about the hit, plus a pointer
// to the actual straw hit.
//
// $Id: StrawHitSynopsis.hh,v 1.1 2011/06/05 23:11:35 mf Exp $
// $Author: mf $
// $Date: 2011/06/05 23:11:35 $
//
// Original author: Mark Fischler
//
//=============================================================================

// C++ includes
#include <iostream>
#include <vector>

// Framework includes.

// Mu2e includes.

namespace mu2e {

  class StrawHitSynopsis 
  {
  public:
    // ctor -- provide the TTracker geometry, the hit, and time estimate
    StrawHitSynopsis(StrawHit const& h, double t0_star)
      : _strawHitPtr (&h)
      , _t0_star (t0_star) {}
      
    StrawHitSynopsis() 
    : _strawHitPtr (0) 
    , _t0_star (0) {}
   
    // TODO - There will be more content to this class; also more accessors
 
    // accessors
    
    double t0_star() const {return _t0_star;}
    StrawHit const * strawHitPtr() const {return _strawHitPtr;}
        
  private:
   // private data
    StrawHit const * _strawHitPtr;
    double _t0_star;
    
  }; // StrawHitSynopsis class definition

  inline
  std::ostream & 
  operator << ( std::ostream & os, const StrawHitSynopsis & s )
  {
    os << "StrawHitSynopsis: straw " <<  s.strawHitPtr()->strawIndex().asInt()
       << " at t0_star = " << s.t0_star() << "\n";
    return os;
  }

} // end of namespace mu2e

#endif // EarlyPatRecPR_StrawHitSynopsis_hh
