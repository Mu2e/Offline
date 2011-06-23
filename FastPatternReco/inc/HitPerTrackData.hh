#ifndef HitPerTrackData_HH
#define HitPerTrackData_HH
//
// this is a old version of Data of the Electrons tracks that came from the targets
//
// $Id: HitPerTrackData.hh,v 1.1 2011/06/23 21:56:11 tassiell Exp $
// $Author: tassiell $
// $Date: 2011/06/23 21:56:11 $
//
// Original author G. Tassielli
//

// C++ includes.
#include <iostream>
#include <cstddef>

// CLHEP includes.
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {
  struct HitPerTrackData{

    // The actual data for this struct.
    size_t  iHit;
    double  mcHitTime;
    bool    isOverlapped;
    bool    isFirst;
    CLHEP::Hep3Vector hitPoint;
    CLHEP::Hep3Vector hitMomentum;

    HitPerTrackData():
      iHit(0),
      mcHitTime(0.0),
      isOverlapped(false),
      isFirst(false){

    }

    HitPerTrackData( size_t &iHit_, double &mcHitTime_, bool &isOverlapped_, bool &isFirst_, CLHEP::Hep3Vector hitPoint_, CLHEP::Hep3Vector hitMomentum_ ):
            iHit(iHit_),
            mcHitTime(mcHitTime_),
            isOverlapped(isOverlapped_),
            isFirst(isFirst_)
    {
            hitPoint=hitPoint_;
            hitMomentum=hitMomentum_;
    }

    // Compiler generated versions are OK for:
    // destructor, copy c'tor, assignment operator.

  };

  inline bool operator==(const HitPerTrackData& lhs,
                         const HitPerTrackData& rhs){
      return ( lhs.iHit == rhs.iHit && lhs.mcHitTime == rhs.mcHitTime &&
               lhs.isOverlapped == rhs.isOverlapped && lhs.isFirst == rhs.isFirst &&
               lhs.hitPoint == rhs.hitPoint && lhs.hitMomentum == rhs.hitMomentum );
  }

  inline bool operator!=(const HitPerTrackData& lhs,
                         const HitPerTrackData& rhs){
    return !(lhs==rhs);
  }

  // Sort first on ProductID and then on index.
  inline bool operator<(const HitPerTrackData& lhs,
                        const HitPerTrackData& rhs){
    return ( lhs.mcHitTime < rhs.mcHitTime );
  }

  // ProductID does not define operators >, <=, >= so we would need
  // to fix that before defining those operators for this class.


  inline std::ostream& operator<<( std::ostream& ost,
                                   HitPerTrackData const& dpi ){
    ost << "i-th hit in the hit sequence " << dpi.iHit
        << " mc hit time "  << dpi.mcHitTime
        << " more tracks share the same hit? "<<dpi.isOverlapped
        << " this track is producing the measure? "<<dpi.isFirst
        <<std::endl
        << " hit position "<<dpi.hitPoint
        << " hit momentum "<<dpi.hitMomentum
        <<std::endl;
    return ost;
  }


}

#endif
