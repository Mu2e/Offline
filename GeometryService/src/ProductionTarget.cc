#include "GeometryService/inc/ProductionTarget.hh"

namespace mu2e {
  ProductionTarget::ProductionTarget(double rOut, double halfLength, double rotX, double rotY, const CLHEP::Hep3Vector& position)
    : _protonBeamRotation(CLHEP::HepRotation::IDENTITY)
    , _prodTargetPosition(position)
    , _rOut(rOut)
    , _halfLength(halfLength)
  {
    _protonBeamRotation.rotateX(rotX).rotateY(rotY);
    _protonBeamInverseRotation = _protonBeamRotation.inverse();
  }
}
