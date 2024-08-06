// Magnet parameters, used by both filter and detector magnets.
//
// Andrei Gaponenko, 2011

#ifndef EXTMONFNALMAGNET_HH
#define EXTMONFNALMAGNET_HH

#include <vector>

#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

namespace mu2e {

  class ExtMonFNALMagnetMaker;

  class ExtMonFNALMagnet {
    std::vector<double> outerHalfSize_;
    double apertureWidth_;
    double apertureHeight_;

    CLHEP::Hep3Vector bfield_;
    double magneticLength_;

    double nominalMomentum_;

    CLHEP::Hep3Vector  refPointInMu2e_;
    CLHEP::HepRotation inRotationInMu2e_;
    CLHEP::HepRotation outRotationInMu2e_;
    CLHEP::HepRotation magnetRotationInMu2e_;

    CLHEP::Hep3Vector geometricCenterInMu2e_;

  public:

    ExtMonFNALMagnet();
    // An initialized instance of this class should be obtained via ExtMonFNALMagnetMaker
    friend class ExtMonFNALMagnetMaker;

    const std::vector<double>& outerHalfSize() const { return outerHalfSize_; }
    double apertureWidth() const { return apertureWidth_; }
    double apertureHeight() const { return apertureHeight_; }

    const CLHEP::Hep3Vector& bfield() const { return bfield_; }
    double magneticLength() const { return magneticLength_; }

    double nominalMomentum() const { return nominalMomentum_; }

    // derived:
    double trackBendRadius(double momentum) const;
    double trackBendHalfAngle(double momentum) const;
    double trackPinvFromRinv(double rinv) const { return trackBendRadius(rinv); }

    double nominalBendHalfAngle() const { return trackBendHalfAngle(nominalMomentum()); }

    // placement
    const CLHEP::Hep3Vector&  refPointInMu2e() const { return refPointInMu2e_; }
    const CLHEP::HepRotation& inRotationInMu2e()     const { return inRotationInMu2e_; }
    const CLHEP::HepRotation& outRotationInMu2e()    const { return outRotationInMu2e_; }
    const CLHEP::HepRotation& magnetRotationInMu2e() const { return magnetRotationInMu2e_; }

    const CLHEP::Hep3Vector&  geometricCenterInMu2e() const { return geometricCenterInMu2e_; }

    // another way to look at rotation
    CLHEP::Hep2Vector dxdzdydz() const;
  };

}// namespace mu2e

#endif/*EXTMONFNALMAGNET_HH*/
