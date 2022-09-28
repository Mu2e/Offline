#include "Offline/PTMGeom/inc/PTM.hh"

// ProductionTarget Monitor(PTM) Object
//
// Author: Helenka Casler
//

namespace mu2e{
  // version 2, includes stand
  PTM::PTM(int version,
          CLHEP::Hep3Vector const& originInMu2e,
          CLHEP::HepRotation const& rotationInMu2e,
          std::shared_ptr<PTMStand> ptmStand,
          std::shared_ptr<PTMHead> ptmHead,
          double totalLength,
          double totalWidth,
          double totalHeight)
          : _version(version),
            _originInMu2e(originInMu2e),
            _rotationInMu2e(rotationInMu2e),
            _ptmStand(ptmStand),
            _ptmHead(ptmHead),
            _totalLength(totalLength),
            _totalWidth(totalWidth),
            _totalHeight(totalHeight) {}

  // version 1, no stand
  PTM::PTM(int version,
          CLHEP::Hep3Vector const& originInMu2e,
          CLHEP::HepRotation const& rotationInMu2e,
          std::shared_ptr<PTMPWC> nearPWC,
          std::shared_ptr<PTMPWC> farPWC,
          double pwcSeparation,
          double motherMargin)
          : _version(version),
            _originInMu2e(originInMu2e),
            _rotationInMu2e(rotationInMu2e),
            _nearPWC(nearPWC),
            _farPWC(farPWC) {
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
