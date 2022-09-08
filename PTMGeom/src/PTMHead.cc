#include "Offline/PTMGeom/inc/PTMHead.hh"

// ProductionTarget Monitor Head (PTM) Object
//
// Author: Helenka Casler
//

namespace mu2e{
  PTMHead::PTMHead(CLHEP::Hep3Vector const& originInMu2e,
          CLHEP::HepRotation const& rotationInMu2e,
          std::shared_ptr<PTMPWC> nearPWC,
          std::shared_ptr<PTMPWC> farPWC,
          double pwcSeparation,
          std::shared_ptr<Box> holderExtrusionLong,
          std::shared_ptr<Box> holderExtrusionShort,
          std::string holderExtrusionMaterialName,
          double holderExtrusionLongSep,
          double holderExtrusionShortPos,
          std::shared_ptr<Box> handleBase,
          std::string handleMaterialName,
          double handleHoleSemiMajor,
          double handleHoleSemiMinor,
          CLHEP::Hep3Vector handleHoleCenter,
          double handleCornerCutSide,
          double motherMargin,
          CLHEP::Hep3Vector const& parentOriginInMu2e,
          CLHEP::HepRotation const& parentRotationInMu2e)
          : _originInMu2e(originInMu2e),
            _rotationInMu2e(rotationInMu2e),
            _nearPWC(nearPWC),
            _farPWC(farPWC),
            _holderExtrusionLong(holderExtrusionLong),
            _holderExtrusionShort(holderExtrusionShort),
            _holderExtrusionMaterialName(holderExtrusionMaterialName),
            _holderExtrusionLongSep(holderExtrusionLongSep),
            _holderExtrusionShortPos(holderExtrusionShortPos),
            _handleBase(handleBase),
            _handleMaterialName(handleMaterialName),
            _handleHoleSemiMajor(handleHoleSemiMajor),
            _handleHoleSemiMinor(handleHoleSemiMinor),
            _handleHoleCenter(handleHoleCenter),
            _handleCornerCutSide(handleCornerCutSide) {

    _totalHeight = holderExtrusionLongSep + (2*handleBase->getYhalfLength()) + motherMargin;
    _totalWidth = holderExtrusionLongSep + (2*holderExtrusionLong->getYhalfLength()) + motherMargin;

    _totalLength = pwcSeparation + 0.5*nearPWC->totalThick() + 0.5*farPWC->totalThick() + motherMargin;

    _originInParent = originInMu2e - parentOriginInMu2e;
    _rotationInParent = rotationInMu2e * parentRotationInMu2e.inverse();

  }

}
