// Andrei Gaponenko, 2011

#ifndef EXTMONFNALMAGNETMAKER_HH
#define EXTMONFNALMAGNETMAKER_HH

#include <string>

#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALMagnet.hh"

namespace CLHEP { class Hep3Vector; }
namespace CLHEP { class HepRotation; }

namespace mu2e {
  class SimpleConfig;

  class ExtMonFNALMagnetMaker {
  public:

    // Initialize sizes and B field, but not orientation and position
    static ExtMonFNALMagnet readIntrinsicParameters(const SimpleConfig& c,
                                                    const std::string& prefix);
    // Initialize position and rotation
    static void positionMagnetRelative(ExtMonFNALMagnet& magnet,
                                       const CLHEP::HepRotation& magnetInRotationInMu2e, // of the input arm of ref trajectory
                                       const CLHEP::Hep3Vector& refTrajMagnetEntranceInMu2e,
                                       double nominalMomentum);

    // This combines the two functions above
    static ExtMonFNALMagnet read(const SimpleConfig& c,
                                 const std::string& prefix,
                                 const CLHEP::HepRotation& magnetInRotationInMu2e, // of the input arm of ref trajectory
                                 const CLHEP::Hep3Vector& refTrajMagnetEntranceInMu2e,
                                 double nominalMomentum);
  };
}

#endif/*EXTMONFNALMAGNETMAKER_HH*/
