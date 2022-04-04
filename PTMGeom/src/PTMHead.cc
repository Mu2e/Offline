#include "Offline/PTMGeom/inc/PTMHead.hh"

// ProductionTarget Monitor Head (PTM) Object
//
// Author: Helenka Casler
//

namespace mu2e{
  PTM::PTMHead(CLHEP::Hep3Vector const& originInMu2e, 
          CLHEP::HepRotation const& rotationInMu2e, 
          std::shared_ptr<PTMPWC> nearPWC, 
          std::shared_ptr<PTMPWC> farPWC,
          double pwcSeparation,
          std::shared_ptr<Box> holderExtrusionLong,
          std::shared_ptr<Box> holderExtrusionShort,
          std::string holderExtrusionMaterialName,
          double holderExtrusionLongSep,
          double holderExtrusionShortPos,
          double motherMargin)
          : _originInMu2e(originInMu2e),
            _rotationInMu2e(rotationInMu2e),
            _nearPWC(nearPWC),
            _farPWC(farPWC),
            _holderExtrusionLong(holderExtrusionLong),
            _holderExtrusionShort(holderExtrusionShort),
            _holderExtrusionMaterialName(holderExtrusionMaterialName),
            _holderExtrusionLongSep(holderExtrusionLongSep),
            _holderExtrusionShortPos(holderExtrusionShortPos) {
    // set the "overall" dimensions based on the PWC dimensions and positions.
    // Not assuming the two PWCs are the same size.
    if (nearPWC->totalHeight() >= farPWC->totalHeight()) {
      _totalHeight = nearPWC->totalHeight()+motherMargin;
    } else {
      _totalHeight = farPWC->totalHeight()+motherMargin;
    }
    if (nearPWC->totalWidth() >= farPWC->totalWidth()) {
      _totalWidth = nearPWC->totalWidth()+motherMargin;
    } else {
      _totalWidth = farPWC->totalWidth()+motherMargin;
    }
    _totalLength = pwcSeparation + 0.5*nearPWC->totalThick() + 0.5*farPWC->totalThick() + motherMargin;

  }

}
