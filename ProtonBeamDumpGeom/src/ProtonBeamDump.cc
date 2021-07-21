#include "Offline/ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"

namespace mu2e {

  ProtonBeamDump::ProtonBeamDump()
    : _coreRotY(0.)
    , _coreRotationInMu2e(CLHEP::HepRotation::IDENTITY)
    , _minCoreShieldingThickness(0.)
    , _shieldingFaceXmin(0.)
    , _shieldingFaceXmax(0.)
    , _shieldingFaceZatXmin(0.)
    , _shieldingFaceZatXmax(0.)
  {}

  //================================================================
  CLHEP::Hep3Vector ProtonBeamDump::mu2eToBeamDump_position(const CLHEP::Hep3Vector& mu2epos) const {

    const CLHEP::Hep3Vector rel(mu2epos - _coreCenterInMu2e);

    static const CLHEP::HepRotation invRot(_coreRotationInMu2e.inverse());
    const CLHEP::Hep3Vector res = invRot * rel;

    // AGDEBUG("POS: mu2e = "<<mu2epos<<", rel = "<<rel<<", res = "<<res);

    return res;
  }

  //================================================================
  CLHEP::Hep3Vector ProtonBeamDump::mu2eToBeamDump_momentum(const CLHEP::Hep3Vector& mu2emom) const {
    static const CLHEP::HepRotation invRot(_coreRotationInMu2e.inverse());
    const CLHEP::Hep3Vector res = invRot * mu2emom;
    // AGDEBUG("MOM: mu2e = "<<mu2emom<<", res = "<<res);
    return res;
  }

  //================================================================
  CLHEP::Hep3Vector ProtonBeamDump::beamDumpToMu2e_position(const CLHEP::Hep3Vector& dumppos) const {
    return _coreCenterInMu2e + _coreRotationInMu2e * dumppos;
  }

  //================================================================
  CLHEP::Hep3Vector ProtonBeamDump::beamDumpToMu2e_momentum(const CLHEP::Hep3Vector& dumpmom) const {
    return _coreRotationInMu2e * dumpmom;
  }

  //================================================================
}
