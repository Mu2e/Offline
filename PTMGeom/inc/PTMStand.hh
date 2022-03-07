#ifndef PTMGeom_PTMStand_hh
#define PTMGeom_PTMStand_hh

#include <memory>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

// ProductionTarget Monitor Stand (PTM) 
//
// Author: Helenka Casler
//

// TODO: ADD THE FOLLOWING:
// - aluminum base that splays out from the center column
// - center column that keeps it at the right height
// - support blocks on the center column
// - top wedge that keeps the detectors at the right angle

namespace mu2e {

  class PTMStand {

  public:
    PTMStand(CLHEP::Hep3Vector const& originInMu2e, 
          CLHEP::HepRotation const& rotationInMu2e);
    PTMStand() {}

    CLHEP::Hep3Vector const &  originInMu2e()   const { return _originInMu2e; }
    CLHEP::HepRotation const & rotationInMu2e() const { return _rotationInMu2e; }




  private:
    CLHEP::Hep3Vector _originInMu2e;
    CLHEP::HepRotation _rotationInMu2e;



  }; // class PTMStand

} // namespace mu2e

#endif