#include "ExtinctionMonitorFNAL/inc/ExtMonFNALBuilding.hh"

#include <cmath>

#include "cetlib/exception.h"

namespace mu2e {
  ExtMonFNALBuilding::ExtMonFNALBuilding()
    : _filterEntranceOffsetX(0.)
    , _filterEntranceOffsetY(0.)
    , _filterAngleH(0.)
    , _filterEntranceAngleV(0.)
    , _extMonFNAL_nominalMomentum(0.)
    , roomInsideFullHeight_(0.)
    , roomWallThickness_(0.)
    , roomFloorThickness_(0.)
    , roomCeilingThickness_(0.)
    , magnetRoomLength_(0.)
    , coll2ShieldingDumpXmin_(0.)
    , coll2ShieldingDumpXmax_(0.)
    , roomInsideYmin_(0.)
    , roomInsideYmax_(0.)
  {}

  //================================================================
  double ExtMonFNALBuilding::FilterMagnetExtMonFNAL::trackBendHalfAngle(double momentum) const {

    // In the bend plane: compute the gyroradius
    // The constant factor is 1/c_light scaled such as
    // to get rTrack in millimeters
    const double rTrack = 3335.64095198 * (momentum/CLHEP::GeV) / (_fieldStrength/CLHEP::tesla);

    //    std::cerr<<"AG: got rTrack = "<<rTrack<<" mm for p = "
    //           <<(momentum/CLHEP::GeV)<<" GeV and  B = "
    //           <<(_fieldStrength()/CLHEP::tesla)<<" tesla"<<std::endl;

    // Can't do momenta that are too low.  For simplicity we just
    // check for the "absolutely impossible" requests here.  The real
    // momentum constraint is tighter because of other pieces of
    // geometry.

    if(_outerHalfSize[2] < rTrack) {
      return asin(_outerHalfSize[2]/rTrack);
    }
    else {
      throw cet::exception("GEOM")<<"ExtMonFNALBuilding::FilterMagnetExtMonFNAL::trackBendHalfAngle(): "
                                  <<"requested momentum p="<<momentum/CLHEP::GeV<<" GeV is too low ";
    }
  }

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
