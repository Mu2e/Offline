// Andrei Gaponenko, 2011

#include "Offline/GeometryService/inc/ExtMonFNALMagnetMaker.hh"

#include <algorithm>
#include <iterator>
#include <iostream>
#include <cmath>

#include "cetlib_except/exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "Offline/ConfigTools/inc/SimpleConfig.hh"

//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

namespace mu2e {

  //================================================================
  ExtMonFNALMagnet
  ExtMonFNALMagnetMaker::readIntrinsicParameters(const SimpleConfig& c,
                                                 const std::string& prefix)
  {
    ExtMonFNALMagnet mag;

    c.getVectorDouble(prefix + ".outerHalfSize", mag.outerHalfSize_, 3);
    mag.apertureWidth_ = c.getDouble(prefix + ".apertureWidth") * CLHEP::mm;
    mag.apertureHeight_ = c.getDouble(prefix + ".apertureHeight") * CLHEP::mm;
    mag.fieldIntegral_ = c.getDouble(prefix + ".BdL") * CLHEP::tesla * CLHEP::mm;
    mag.magneticLength_ = c.getDouble(prefix + ".magneticLength") * CLHEP::mm;

    return mag;
  }

  //================================================================
  void ExtMonFNALMagnetMaker::
  positionMagnetRelative(ExtMonFNALMagnet& mag,
                         const CLHEP::HepRotation& magnetInRotationInMu2e, // of the input arm of ref trajectory
                         const CLHEP::Hep3Vector& refTrajMagnetEntranceInMu2e,
                         double nominalMomentum)
  {
    mag.nominalMomentum_ = nominalMomentum;
    mag.inRotationInMu2e_ = magnetInRotationInMu2e;

    const double fieldStrength = mag.fieldIntegral_ / mag.magneticLength_;
    mag.bfield_ = mag.inRotationInMu2e_ * CLHEP::Hep3Vector(fieldStrength, 0, 0);

    const double trackBendHalfAngle = mag.trackBendHalfAngle(nominalMomentum);
    mag.outRotationInMu2e_ = CLHEP::HepRotation(mag.bfield_, -2*trackBendHalfAngle) * magnetInRotationInMu2e;
    mag.magnetRotationInMu2e_ = CLHEP::HepRotation(mag.bfield_, -trackBendHalfAngle) * magnetInRotationInMu2e;

    // distance between the points
    const double magnetEntranceToBendPointDistance = mag.outerHalfSize_[2]/cos(trackBendHalfAngle);

    mag.refPointInMu2e_ = refTrajMagnetEntranceInMu2e
      + magnetInRotationInMu2e * CLHEP::Hep3Vector(0,0, -magnetEntranceToBendPointDistance);

    // Putting the geometric center at this distance from the bend point
    // would center the "kinked" trajectory in the magnet
    const double refCenterDistanceKinked = mag.outerHalfSize_[2] * tan(trackBendHalfAngle)/2;

    // For magneticLength != 0 the trajectory has a circular part instead ofthe kink
    // Distance between the reference point and the circle
    const double refCenterDistanceToTrajectory = (mag.magneticLength_/2)*(1./cos(trackBendHalfAngle) - 1)/sin(trackBendHalfAngle);

    const double bendToCenterDistance = refCenterDistanceKinked + 0.5 * refCenterDistanceToTrajectory;

    mag.geometricCenterInMu2e_ = mag.refPointInMu2e_ + mag.magnetRotationInMu2e_ * CLHEP::Hep3Vector(0, -bendToCenterDistance, 0);
  }

  //================================================================
  ExtMonFNALMagnet
  ExtMonFNALMagnetMaker::read(const SimpleConfig& c, const std::string& prefix,
                              const CLHEP::HepRotation& magnetInRotationInMu2e, // of the input arm of ref trajectory
                              const CLHEP::Hep3Vector& refTrajMagnetEntranceInMu2e,
                              double nominalMomentum)
  {
    ExtMonFNALMagnet mag = readIntrinsicParameters(c, prefix);
    positionMagnetRelative(mag, magnetInRotationInMu2e, refTrajMagnetEntranceInMu2e, nominalMomentum);

    if(c.getInt("extMonFNAL.verbosityLevel") > 0) {
      std::cout<<"ExtMonFNALMagnet "<<prefix<<": refPointInMu2e = "<<mag.refPointInMu2e_<<std::endl;
      std::cout<<"ExtMonFNALMagnet "<<prefix<<": geometricCenterInMu2e = "<<mag.geometricCenterInMu2e_<<std::endl;
    }

    return mag;
  }

  //================================================================

} // namespace mu2e
