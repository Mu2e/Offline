#include "ProductionTargetGeom/inc/ProductionTarget.hh"

namespace mu2e {
  ProductionTarget::ProductionTarget(int version, double rOut, 
				     double halfLength, double rotX, 
				     double rotY, const CLHEP::Hep3Vector& position,
				     int    nFins = 0,
				     double finHt = 0, double finThick = 0,
				     double hubDisU = 0, double hubDisD = 0,
				     double hubAngU = 0, double hubAngD = 0,
				     double huboverU = 0, double huboverD = 0)
    : _protonBeamRotation(CLHEP::HepRotation::IDENTITY)
    , _prodTargetPosition(position)
    , _version(version)
    , _rOut(rOut)
    , _halfLength(halfLength)
    , _envelHalfLength(halfLength)
    , _finHeight(finHt)
    , _finThickness(finThick)
    , _hubDistUS(hubDisU)
    , _hubDistDS(hubDisD)
    , _hubAngleUS(hubAngU)
    , _hubAngleDS(hubAngD)
    , _hubOverhangUS(huboverU)
    , _hubOverhangDS(huboverD)
  {
    _protonBeamRotation.rotateX(rotX).rotateY(rotY);
    _protonBeamInverseRotation = _protonBeamRotation.inverse();
  }
}
