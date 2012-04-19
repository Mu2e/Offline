// Andrei Gaponenko, 2011

#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL.hh"

#include <iostream>

//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    ExtMon::ExtMon(const std::vector<double>& logicalEnclosureHalfDim,
                   double roomCenterHeightAboveDumpCore)

      : m_roomHalfSize(logicalEnclosureHalfDim)
      , m_roomCenterHeightAboveDumpCore(roomCenterHeightAboveDumpCore)
    {}

    //================================================================
    std::vector<double> ExtMon::sensorHalfSize(unsigned iplane) const {
      std::vector<double> res(3);
      res[0] = sensor_halfdx()[iplane];
      res[1] = sensor_halfdy()[iplane];
      res[2] = sensor_halfdz()[iplane];
      return res;
    }

    //================================================================
    CLHEP::Hep3Vector ExtMon::sensorOffsetInParent(unsigned iplane) const {
      CLHEP::Hep3Vector res(3);
      res[0] = sensor_xoffset()[iplane];
      res[1] = sensor_yoffset()[iplane];
      res[2] = sensor_zoffset()[iplane];
      return res;
    }

    //================================================================
    CLHEP::Hep3Vector ExtMon::mu2eToExtMon_position(const CLHEP::Hep3Vector& mu2epos) const {
      const CLHEP::Hep3Vector rel(mu2epos - m_coordinateCenterInMu2e);
      const CLHEP::Hep3Vector res = m_coordinateRotationInMu2e * rel;
      AGDEBUG("POS: mu2e = "<<mu2epos<<", rel = "<<rel<<", res = "<<res);
      return res;
    }

    //================================================================
    CLHEP::Hep3Vector ExtMon::mu2eToExtMon_momentum(const CLHEP::Hep3Vector& mu2emom) const {
      const CLHEP::Hep3Vector res = m_coordinateRotationInMu2e * mu2emom;
      AGDEBUG("MOM: mu2e = "<<mu2emom<<", res = "<<res);
      return res;
    }

    //================================================================

  } // namespace ExtMonFNAL
} // namespace mu2e
