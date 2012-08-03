// Andrei Gaponenko, 2011

#ifndef EXTMONFNALMAGNETMAKER_HH
#define EXTMONFNALMAGNETMAKER_HH

#include <string>

#include "ExtinctionMonitorFNAL/inc/ExtMonFNALMagnet.hh"

namespace CLHEP { class Hep3Vector; }
namespace CLHEP { class HepRotation; }

namespace mu2e {
  class SimpleConfig;

  class ExtMonFNALMagnetMaker {
  public:

    static ExtMonFNALMagnet read(const SimpleConfig& c,
                                 const std::string& prefix,
                                 const CLHEP::Hep3Vector& magnetRefPointInMu2e,
                                 const CLHEP::HepRotation& magnetInRotation,
                                 double nominalMomentum);
  };
}

#endif/*EXTMONFNALMAGNETMAKER_HH*/
