#include "PTMonGeom/inc/PTMon.hh"

// ProductionTarget Monitor (PTMon) Object
//
// Author: Helenka Casler
//

namespace mu2e{
  PTMon::PTMon(CLHEP::Hep3Vector const& originInMu2e, 
          CLHEP::HepRotation const& rotationInMu2e, 
          PTMonPWC const& nearPWC, PTMonPWC const& farPWC,
          double pwcSeparation)
          : _originInMu2e(originInMu2e),
            _rotationInMu2e(rotationInMu2e),
            _nearPWC(nearPWC),
            _farPWC(farPWC) {
    // set the "overall" dimensions based on the PWC dimensions and positions.
    // Not assuming the two PWCs are the same size.
    if (nearPWC->frameHeight() >= farPWC->frameHeight()) {
      _totalHeight = nearPWC->frameHeight();
    } else {
      _totalHeight = farPWC->frameHeight();
    }
    if (nearPWC->frameWidth() >= farPWC->frameWidth()) {
      _totalWidth = nearPWC->frameWidth();
    } else {
      _totalWidth = farPWC->frameWidth();
    }
    _totalLength = pwcSeparation + 0.5*nearPWC->totalThick() + 0.5*farPWC->totalThick();

  }

}