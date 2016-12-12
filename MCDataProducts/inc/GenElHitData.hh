#ifndef MCDataProducts_GenElHitData_hh
#define MCDataProducts_GenElHitData_hh
//
// all the hits data of the Electrons tracks that came from the targets
//
// $Id: GenElHitData.hh,v 1.1 2011/10/11 17:27:52 tassiell Exp $
// $Author: tassiell $
// $Date: 2011/10/11 17:27:52 $
//
// Original author G. Tassielli
//


// C++ includes.
#include <iostream>
#include <cstddef>

// CLHEP includes.
#include "CLHEP/Vector/ThreeVector.h"

// Mu2e includes.
#include "RecoDataProducts/inc/StrawHit.hh"
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  typedef art::Ptr<StrawHit> StrawHitPtr;

  struct GenElHitData{

    // The actual data for this struct.
    StrawHitPtr  _iHit;
    double       _mcHitTime;
    bool         _isOverlapped;
    bool         _isOvrlpByProton;
    bool         _isFirst;
    CLHEP::Hep3Vector _hitPoint;
    CLHEP::Hep3Vector _hitMomentum;

    GenElHitData():
      _mcHitTime(0.00000),
      _isOverlapped(false),
      _isOvrlpByProton(false),
      _isFirst(false){

    }

    GenElHitData( StrawHitPtr &iHit_, double &mcHitTime_, bool &isOverlapped_, bool &isOvrlpByProton_, bool &isFirst_, CLHEP::Hep3Vector hitPoint_, CLHEP::Hep3Vector hitMomentum_ ):
            _iHit(iHit_),
            _mcHitTime(mcHitTime_),
            _isOverlapped(isOverlapped_),
            _isOvrlpByProton(isOvrlpByProton_),
            _isFirst(isFirst_)
    {
            _hitPoint=hitPoint_;
            _hitMomentum=hitMomentum_;
    }

    // Compiler generated versions are OK for:
    // destructor, copy c'tor, assignment operator.

  };

  inline bool operator==(const GenElHitData& lhs,
                         const GenElHitData& rhs){
      return ( lhs._iHit.key() == rhs._iHit.key() && lhs._mcHitTime == rhs._mcHitTime &&
               lhs._isOverlapped == rhs._isOverlapped && lhs._isOvrlpByProton == rhs._isOvrlpByProton && lhs._isFirst == rhs._isFirst &&
               lhs._hitPoint == rhs._hitPoint && lhs._hitMomentum == rhs._hitMomentum );
  }

  inline bool operator!=(const GenElHitData& lhs,
                         const GenElHitData& rhs){
    return !(lhs==rhs);
  }

  // Sort first on ProductID and then on index.
  inline bool operator<(const GenElHitData& lhs,
                        const GenElHitData& rhs){
    return ( lhs._mcHitTime < rhs._mcHitTime );
  }

  // ProductID does not define operators >, <=, >= so we would need
  // to fix that before defining those operators for this class.


  inline std::ostream& operator<<( std::ostream& ost,
                                   GenElHitData const& dpi ){
    ost << "i-th hit in the hit sequence " << dpi._iHit.key()
        << " mc hit time "  << dpi._mcHitTime
        << " more tracks share the same hit? "<<dpi._isOverlapped<<" if yes is it a proton? "<<dpi._isOvrlpByProton
        << " this track is producing the measure? "<<dpi._isFirst
        <<std::endl
        << " hit position "<<dpi._hitPoint
        << " hit momentum "<<dpi._hitMomentum
        <<std::endl;
    return ost;
  }


} // namespace mu2e

#endif /* MCDataProducts_GenElHitData_hh */
