// Andrei Gaponenko, 2011

#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL.hh"

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    ExtMon::ExtMon(const std::vector<double>& logicalEnclosureHalfDim,
                   double roomCenterHeightAboveDumpCore)

      : m_roomHalfSize(logicalEnclosureHalfDim)
      , m_roomCenterHeightAboveDumpCore(roomCenterHeightAboveDumpCore)
      , m_detectorRotationInRoom(CLHEP::HepRotation::IDENTITY)

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

  }
}
