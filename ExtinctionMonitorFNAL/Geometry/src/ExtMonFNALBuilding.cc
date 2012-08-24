#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"

#include <cmath>

#include "cetlib/exception.h"

namespace mu2e {
  ExtMonFNALBuilding::ExtMonFNALBuilding()
    : _filterEntranceOffsetX(0.)
    , _filterEntranceOffsetY(0.)
    , _filterAngleH(0.)
    , _filterEntranceAngleV(0.)
    , roomInsideFullHeight_(0.)
    , roomWallThickness_(0.)
    , roomFloorThickness_(0.)
    , roomCeilingThickness_(0.)
    , dirtOverheadThickness_(0.)
    , dirtOverheadHorizontalMargin_(0.)
    , magnetRoomLength_(0.)
    , coll2ShieldingDumpXmin_(0.)
    , coll2ShieldingDumpXmax_(0.)
    , roomInsideYmin_(0.)
    , roomInsideYmax_(0.)
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
  CLHEP::Hep3Vector ExtMonFNALBuilding::ceilingRefPointInMu2e() const {
    return roomRefPointInMu2e() + CLHEP::Hep3Vector
      (0, 0.5*(roomInsideFullHeight_ + roomCeilingThickness_), 0);
  }

  //================================================================
  CLHEP::Hep3Vector ExtMonFNALBuilding::floorRefPointInMu2e() const {
    return roomRefPointInMu2e() + CLHEP::Hep3Vector
      (0, -0.5*(roomInsideFullHeight_ + roomFloorThickness_), 0);
  }

  //================================================================

}
