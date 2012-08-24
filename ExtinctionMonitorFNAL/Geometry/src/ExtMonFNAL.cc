// Andrei Gaponenko, 2011

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"

#include <iostream>

//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    CLHEP::Hep3Vector ExtMon::detectorCenterInMu2e() const {
      return up_.refPointInMu2e();
    }

    //================================================================
    const CLHEP::HepRotation& ExtMon::detectorRotationInMu2e() const {
      return up_.rotationInMu2e();
    }

    //================================================================
    CLHEP::Hep3Vector ExtMon::mu2eToExtMon_position(const CLHEP::Hep3Vector& mu2epos) const {
      return up_.mu2eToStack_position(mu2epos);
    }

    //================================================================
    CLHEP::Hep3Vector ExtMon::mu2eToExtMon_momentum(const CLHEP::Hep3Vector& mu2emom) const {
      return up_.mu2eToStack_momentum(mu2emom);
    }

    //================================================================

  } // namespace ExtMonFNAL
} // namespace mu2e
