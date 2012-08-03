// Andrei Gaponenko, 2011

#include "ExtinctionMonitorFNAL/inc/ExtMonFNALMagnetMaker.hh"

#include <algorithm>
#include <iterator>
#include <iostream>
#include <cmath>

#include "cetlib/exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "ConfigTools/inc/SimpleConfig.hh"

//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

namespace mu2e {

  ExtMonFNALMagnet
  ExtMonFNALMagnetMaker::read(const SimpleConfig& c, const std::string& prefix,
                              const CLHEP::Hep3Vector& magnetRefPointInMu2e,
                              const CLHEP::HepRotation& magnetInRotation,
                              double nominalMomentum)
  {
    ExtMonFNALMagnet mag;

    c.getVectorDouble(prefix + ".outerHalfSize", mag.outerHalfSize_, 3);
    mag.apertureWidth_ = c.getDouble(prefix + ".apertureWidth") * CLHEP::mm;
    mag.apertureHeight_ = c.getDouble(prefix + ".apertureHeight") * CLHEP::mm;

    mag.refPointInMu2e_ = magnetRefPointInMu2e;
    mag.inRotationInMu2e_ = magnetInRotation;
    mag.nominalMomentum_ = nominalMomentum;

    const double fieldStrength = c.getDouble(prefix + ".fieldStrength") * CLHEP::tesla;
    mag.bfield_ = mag.inRotationInMu2e_ * CLHEP::Hep3Vector(fieldStrength, 0, 0);

    const double halfAlpha = mag.trackBendHalfAngle(nominalMomentum);

    mag.outRotationInMu2e_ = CLHEP::HepRotation(mag.bfield_, -2*halfAlpha) * magnetInRotation;
    mag.magnetRotationInMu2e_ = CLHEP::HepRotation(mag.bfield_, -halfAlpha) * magnetInRotation;

    const double rTrack = mag.trackBendRadius(nominalMomentum);
    const double trackBendHalfAngle = mag.trackBendHalfAngle(nominalMomentum);
    const double refCenterDistance = rTrack*(1/cos(trackBendHalfAngle) - 0.5*(1+cos(trackBendHalfAngle)));

    mag.geometricCenterInMu2e_ = mag.refPointInMu2e_ + mag.magnetRotationInMu2e_ * CLHEP::Hep3Vector(0, -refCenterDistance, 0);

    AGDEBUG("refCenterDistance = "<<refCenterDistance);
    AGDEBUG("refPointInMu2e = "<<mag.refPointInMu2e_);
    AGDEBUG("geometricCenterInMu2e = "<<mag.geometricCenterInMu2e_);

    return mag;
  }

} // namespace mu2e
