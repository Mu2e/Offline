#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"

#include <cmath>

#include "cetlib/exception.h"

namespace mu2e {

  ProtonBeamDump::ProtonBeamDump()
    : _enclosureRotationInMu2e(CLHEP::HepRotation::IDENTITY)
    , _collimator1RotationInMu2e(CLHEP::HepRotation::IDENTITY)
    , _filterMagnetRotationInMu2e(CLHEP::HepRotation::IDENTITY)
    , _collimator2RotationInMu2e(CLHEP::HepRotation::IDENTITY)
  {}

  double ProtonBeamDump::FilterMagnetExtMonFNAL::trackBendHalfAngle(double momentum) const {

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
      throw cet::exception("GEOM")<<"ProtonBeamDump::FilterMagnetExtMonFNAL::trackBendHalfAngle(): "
                                  <<"requested momentum p="<<momentum/CLHEP::GeV<<" GeV is too low ";
    }
  }

  //================================================================
  CLHEP::Hep3Vector ProtonBeamDump::filterEntranceInMu2e() const {
    CLHEP::Hep3Vector filterEntranceInEnclosure(_coreCenterInEnclosure[0] + _filterEntranceOffsetX,
                                                _coreCenterInEnclosure[1] + _filterEntranceOffsetY,
                                                _enclosureHalfSize[2]);

    return _enclosureCenterInMu2e + _enclosureRotationInMu2e * filterEntranceInEnclosure;
  }

  //================================================================
  CLHEP::Hep3Vector ProtonBeamDump::filterExitInMu2e() const {
    CLHEP::Hep3Vector filterExitInEnclosure(_collimator2CenterInEnclosure[0]
                                            + 0.5*_collimator2.horizontalLength()*tan(_collimator2.angleH()),

                                            _collimator2CenterInEnclosure[1]
                                            + 0.5*_collimator2.horizontalLength()*tan(_collimator2.angleV())/cos(_collimator2.angleH()),

                                            -_enclosureHalfSize[2]);

    return _enclosureCenterInMu2e + _enclosureRotationInMu2e * filterExitInEnclosure;
  }

  //================================================================
  CLHEP::Hep3Vector ProtonBeamDump::mu2eToBeamDump_position(const CLHEP::Hep3Vector& mu2epos) const {

    const CLHEP::Hep3Vector rel(mu2epos - _coreCenterInMu2e);

    static const CLHEP::HepRotation invRot(_enclosureRotationInMu2e.inverse());
    const CLHEP::Hep3Vector res = invRot * rel;

    // AGDEBUG("POS: mu2e = "<<mu2epos<<", rel = "<<rel<<", res = "<<res);

    return res;
  }

  //================================================================
  CLHEP::Hep3Vector ProtonBeamDump::mu2eToBeamDump_momentum(const CLHEP::Hep3Vector& mu2emom) const {
    static const CLHEP::HepRotation invRot(_enclosureRotationInMu2e.inverse());
    const CLHEP::Hep3Vector res = invRot * mu2emom;
    // AGDEBUG("MOM: mu2e = "<<mu2emom<<", res = "<<res);
    return res;
  }

  //================================================================
}
