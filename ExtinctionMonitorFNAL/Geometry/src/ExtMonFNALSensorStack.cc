// Andrei Gaponenko, 2012

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALSensorStack.hh"

#include <iostream>

//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

namespace mu2e {

  //================================================================
  ExtMonFNALSensorStack::ExtMonFNALSensorStack()
    : distanceToTestMaterials_()
    , m_testMaterialPitch()
  {}

  //================================================================
  std::vector<double> ExtMonFNALSensorStack::sensorHalfSize(unsigned iplane) const {
    std::vector<double> res(3);
    res[0] = sensor_halfdx()[iplane];
    res[1] = sensor_halfdy()[iplane];
    res[2] = sensor_halfdz()[iplane];
    return res;
  }

  //================================================================
  CLHEP::Hep3Vector ExtMonFNALSensorStack::sensorOffsetInStack(unsigned iplane) const {
    CLHEP::Hep3Vector res(3);
    res[0] = sensor_xoffset()[iplane];
    res[1] = sensor_yoffset()[iplane];
    res[2] = sensor_zoffset()[iplane];
    return res;
  }

  //================================================================
  CLHEP::Hep3Vector ExtMonFNALSensorStack::mu2eToStack_position(const CLHEP::Hep3Vector& mu2epos) const {
    const CLHEP::Hep3Vector rel(mu2epos - m_stackRefPointInMu2e);
    const CLHEP::Hep3Vector res = m_coordinateRotationInMu2e * rel;
    AGDEBUG("POS: mu2e = "<<mu2epos<<", rel = "<<rel<<", res = "<<res);
    return res;
  }

  //================================================================
  CLHEP::Hep3Vector ExtMonFNALSensorStack::mu2eToStack_momentum(const CLHEP::Hep3Vector& mu2emom) const {
    const CLHEP::Hep3Vector res = m_coordinateRotationInMu2e * mu2emom;
    AGDEBUG("MOM: mu2e = "<<mu2emom<<", res = "<<res);
    return res;
  }

  //================================================================

} // namespace mu2e
