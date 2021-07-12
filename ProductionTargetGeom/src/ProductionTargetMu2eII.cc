#include "Offline/ProductionTargetGeom/inc/ProductionTargetMu2eII.hh"

namespace mu2e {
  ProductionTargetMu2eII::ProductionTargetMu2eII(std::string type, int version)
    : _protonBeamRotation(CLHEP::HepRotation::IDENTITY)
    , _prodTargetPosition(CLHEP::Hep3Vector())
    , _type(type)
    , _version(version)
    , _isRotating(false)
    , _isConveyor(false)
  {
    // _protonBeamRotation.rotateX(rotX).rotateY(rotY);
    // _protonBeamInverseRotation = _protonBeamRotation.inverse();
  }

}
