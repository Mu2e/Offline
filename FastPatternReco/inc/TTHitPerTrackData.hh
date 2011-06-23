#ifndef TTHitPerTrackData_HH
#define TTHitPerTrackData_HH
//
// this is a old version of Data of the Electrons tracks that came from the targets
//
// $Id: TTHitPerTrackData.hh,v 1.1 2011/06/23 21:56:11 tassiell Exp $
// $Author: tassiell $
// $Date: 2011/06/23 21:56:11 $
//
// Original author G. Tassielli
//

// Mu2e includes.
#include "FastPatternReco/inc/HitPerTrackData.hh"

namespace mu2e {
  struct TTHitPerTrackData : HitPerTrackData {

    // TT cell id added
    int     strawn;
    int     layern;
    int     sectorn;
    int     devicen;

    TTHitPerTrackData():
      HitPerTrackData(),
      strawn(0),
      layern(0),
      sectorn(0),
      devicen(0){

    }

    TTHitPerTrackData( size_t &iHit_, double &mcHitTime_, bool &isOverlapped_, bool &isFirst_, int &strawn_, int &layern_, int &sectorn_, int &devicen_, CLHEP::Hep3Vector hitPoint_, CLHEP::Hep3Vector hitMomentum_ ):
            HitPerTrackData(iHit_,mcHitTime_,isOverlapped_,isFirst_,hitPoint_,hitMomentum_),
            strawn(strawn_),
            layern(layern_),
            sectorn(sectorn_),
            devicen(devicen_)
    {
    }

    // Compiler generated versions are OK for:
    // destructor, copy c'tor, assignment operator.

  };

  inline bool operator==(const TTHitPerTrackData& lhs,
                         const TTHitPerTrackData& rhs){
      return ( lhs.iHit == rhs.iHit && lhs.mcHitTime == rhs.mcHitTime &&
               lhs.isOverlapped == rhs.isOverlapped && lhs.isFirst == rhs.isFirst &&
               lhs.strawn == rhs.strawn && lhs.layern == rhs.layern && lhs.sectorn == rhs.sectorn && lhs.devicen == rhs.devicen &&
               lhs.hitPoint == rhs.hitPoint && lhs.hitMomentum == rhs.hitMomentum );
  }

  inline bool operator!=(const TTHitPerTrackData& lhs,
                         const TTHitPerTrackData& rhs){
    return !(lhs==rhs);
  }

  // Sort first on ProductID and then on index.
  inline bool operator<(const TTHitPerTrackData& lhs,
                        const TTHitPerTrackData& rhs){
    return ( lhs.mcHitTime < rhs.mcHitTime );
  }

  // ProductID does not define operators >, <=, >= so we would need
  // to fix that before defining those operators for this class.


  inline std::ostream& operator<<( std::ostream& ost,
                                   TTHitPerTrackData const& dpi ){
    ost << "i-th hit in the hit sequence " << dpi.iHit
        << " straw id "<<dpi.sectorn<<" - "<<dpi.layern<<" - "<<dpi.strawn<<" - "<<dpi.devicen
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
