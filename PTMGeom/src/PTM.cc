#include "Offline/PTMGeom/inc/PTM.hh"

// ProductionTarget Monitor(PTM) Object
//
// Author: Helenka Casler
//

// TODO: set overall length, width, height

namespace mu2e{
  PTM::PTM(CLHEP::Hep3Vector const& originInMu2e,
          CLHEP::HepRotation const& rotationInMu2e,
          std::shared_ptr<PTMStand> ptmStand,
          std::shared_ptr<PTMHead> ptmHead)
          : _originInMu2e(originInMu2e),
            _rotationInMu2e(rotationInMu2e),
            _ptmStand(ptmStand),
            _ptmHead(ptmHead) {
    // set the "overall" dimensions based on the dimenions and positions of the
    // head and stand
    _totalHeight = ptmHead->totalHeight() + 4000.; // TODO: make this depend on the stand
    _totalLength = 1220+100;
    _totalWidth = 1220+100;

  }
}
