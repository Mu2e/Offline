// Andrei Gaponenko, 2011

#ifndef EXTMONFNALMAGNETMAKER_HH
#define EXTMONFNALMAGNETMAKER_HH

#include <string>

#include "ExtinctionMonitorFNAL/inc/ExtMonFNALMagnet.hh"

namespace mu2e {
  class SimpleConfig;

  class ExtMonFNALMagnetMaker {
  public:
    static ExtMonFNALMagnet read(const SimpleConfig& c, const std::string& prefix);
  };
}

#endif/*EXTMONFNALMAGNETMAKER_HH*/
