#include "Offline/PTMGeom/inc/PTMStand.hh"

// ProductionTarget Monitor Stand (PTM) Object
//
// Author: Helenka Casler
//

namespace mu2e{
  PTMStand::PTMStand(CLHEP::Hep3Vector const& originInMu2e, 
          CLHEP::HepRotation const& rotationInMu2e)
          : _originInMu2e(originInMu2e),
            _rotationInMu2e(rotationInMu2e) {
    // TODO: fill in the other stuff

  }

}
