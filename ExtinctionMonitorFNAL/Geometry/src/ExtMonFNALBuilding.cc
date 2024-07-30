#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"

#include <cmath>

#include "cetlib_except/exception.h"

namespace mu2e {
  ExtMonFNALBuilding::ExtMonFNALBuilding()
    : _filterAngleH(0.)
    , _filterEntranceAngleV(0.)
    , roomInsideFullHeight_(0.)
    , magnetRoomLength_(0.)
    , roomInsideYmin_(0.)
    , roomInsideYmax_(0.)
    , _HVACductRadius(0.)
    , _HVACductHalfLength(0.)
  {}

  //================================================================
  ExtMonFNALBuilding::CollimatorExtMonFNAL::CollimatorExtMonFNAL()
    : _horizontalLength(0.)
    , _shotLinerOuterRadius(0.)
    , _shotLinerOuterThickness(0.)
    , _radiusTransitiondZ(0.)
    , _angleH(0.)
    , _angleV(0.)
  {}

  //================================================================
  double ExtMonFNALBuilding::CollimatorExtMonFNAL::halfLength() const {
    using std::pow;
    return 0.5*_horizontalLength *
      sqrt(1
           + pow(tan(_angleH), 2)
           + pow(tan(_angleV)/cos(_angleH), 2)
           );
  }

  //================================================================
  CLHEP::Hep3Vector ExtMonFNALBuilding::filterEntranceInMu2e() const {
    return _collimator1CenterInMu2e +
      _collimator1RotationInMu2e * CLHEP::Hep3Vector(0,0, +_collimator1.halfLength());
  }

  //================================================================
  CLHEP::Hep3Vector ExtMonFNALBuilding::filterExitInMu2e() const {
    return _collimator2CenterInMu2e +
      _collimator2RotationInMu2e * CLHEP::Hep3Vector(0,0, -_collimator2.halfLength());
  }

  //================================================================

}
