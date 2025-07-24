#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALCollimator.hh"
#include "Offline/ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"

namespace mu2e {
  //================================================================
  ExtMonFNALCollimator::ExtMonFNALCollimator()
    : _shotLinerOuterRadius(0.)
    , _shotLinerOuterThickness(0.)
    , _radiusTransitiondZ(0.)
    , _angleH_inBeamDump(0.)
    , _angleV(0.)
  {}

  //================================================================
  void ExtMonFNALCollimator::setFromDumpAngles(double angleH_inBeamDump, double angleV, const ProtonBeamDump& dump) {

    _angleH_inBeamDump = angleH_inBeamDump;
    _angleV = angleV;

    CLHEP::HepRotation tmp(CLHEP::HepRotation::IDENTITY);
    tmp.rotateX(angleV).rotateY(-angleH_inBeamDump);
    _rotationInMu2e = dump.coreRotationInMu2e() * tmp;
  }

  //================================================================
  void ExtMonFNALCollimator::setFromMu2eSlopes(double dxdz, double dydz, const ProtonBeamDump& dump) {

    const double angleH_inMu2e = atan(dxdz);
    const double angleV = atan2( -dydz, sqrt(1. + dxdz*dxdz));

    _angleH_inBeamDump = dump.coreRotY() - angleH_inMu2e;
    _angleV = angleV;

    CLHEP::HepRotation tmp(CLHEP::HepRotation::IDENTITY);
    tmp.rotateX(angleV).rotateY(angleH_inMu2e);
    _rotationInMu2e = tmp;
  }

  //================================================================
  CLHEP::Hep3Vector ExtMonFNALCollimator::trajectoryPointInMu2e(double pos) const {
    return _centerInMu2e + _rotationInMu2e * CLHEP::Hep3Vector(0,0, pos);
  }

  //================================================================
  CLHEP::Hep2Vector ExtMonFNALCollimator::dxdzdydz() const {
    const auto v = _rotationInMu2e * CLHEP::Hep3Vector(0,0, -1);
    return CLHEP::Hep2Vector( v.x()/v.z(), v.y()/v.z());
  }

  //================================================================
}
