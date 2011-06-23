#ifndef ITHitPerTrackData_HH
#define ITHitPerTrackData_HH
//
// this is a old version of Data of the Electrons tracks that came from the targets
//
// $Id: ITHitPerTrackData.hh,v 1.1 2011/06/23 21:56:11 tassiell Exp $
// $Author: tassiell $
// $Date: 2011/06/23 21:56:11 $
//
// Original author G. Tassielli
//


// Mu2e includes.
#include "FastPatternReco/inc/HitPerTrackData.hh"

namespace mu2e {
  struct ITHitPerTrackData : HitPerTrackData {

    // IT cell id added
    int     celln;
    int     layern;
    int     slayern;

    ITHitPerTrackData():
      HitPerTrackData(),
      celln(0),
      layern(0),
      slayern(0){

    }

    ITHitPerTrackData( size_t &iHit_, double &mcHitTime_, bool &isOverlapped_, bool &isFirst_, int &celln_, int &layern_, int &slayern_, CLHEP::Hep3Vector hitPoint_, CLHEP::Hep3Vector hitMomentum_ ):
            HitPerTrackData(iHit_,mcHitTime_,isOverlapped_,isFirst_,hitPoint_,hitMomentum_),
            celln(celln_),
            layern(layern_),
            slayern(slayern_)
    {
    }

    // Compiler generated versions are OK for:
    // destructor, copy c'tor, assignment operator.

  };

  inline bool operator==(const ITHitPerTrackData& lhs,
                         const ITHitPerTrackData& rhs){
      return ( lhs.iHit == rhs.iHit && lhs.mcHitTime == rhs.mcHitTime &&
               lhs.isOverlapped == rhs.isOverlapped && lhs.isFirst == rhs.isFirst &&
               lhs.celln == rhs.celln && lhs.layern == rhs.layern && lhs.slayern == rhs.slayern &&
               lhs.hitPoint == rhs.hitPoint && lhs.hitMomentum == rhs.hitMomentum );
  }

  inline bool operator!=(const ITHitPerTrackData& lhs,
                         const ITHitPerTrackData& rhs){
    return !(lhs==rhs);
  }

  // Sort first on ProductID and then on index.
  inline bool operator<(const ITHitPerTrackData& lhs,
                        const ITHitPerTrackData& rhs){
    return ( lhs.mcHitTime < rhs.mcHitTime );
  }

  // ProductID does not define operators >, <=, >= so we would need
  // to fix that before defining those operators for this class.


  inline std::ostream& operator<<( std::ostream& ost,
                                   ITHitPerTrackData const& dpi ){
    ost << "i-th hit in the hit sequence " << dpi.iHit
        << " cell id "<<dpi.slayern<<" - "<<dpi.layern<<" - "<<dpi.celln
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
