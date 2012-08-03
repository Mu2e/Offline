// Andrei Gaponenko, 2011

#include "ExtinctionMonitorFNAL/inc/ExtMonFNALMagnetMaker.hh"

#include <algorithm>
#include <iterator>
#include <iostream>
#include <cmath>

#include "cetlib/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "ConfigTools/inc/SimpleConfig.hh"

#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
//#define AGDEBUG(stuff)

namespace mu2e {

  ExtMonFNALMagnet
  ExtMonFNALMagnetMaker::read(const SimpleConfig& c, const std::string& prefix) {
    ExtMonFNALMagnet mag;
    c.getVectorDouble(prefix + ".outerHalfSize", mag._outerHalfSize, 3);
    mag._apertureWidth = c.getDouble(prefix + ".apertureWidth") * CLHEP::mm;
    mag._apertureHeight = c.getDouble("extMonFNAL.magnet.apertureHeight") * CLHEP::mm;
    mag._fieldStrength = c.getDouble("extMonFNAL.magnet.fieldStrength") * CLHEP::tesla;
    return mag;
  }

} // namespace mu2e
