// Sam Fine, 2024

#ifndef CONSTRUCTEXTMONFNALInfrastructure_HH
#define CONSTRUCTEXTMONFNALInfrastructure_HH

#include <string>

#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPlaneStack.hh"

namespace CLHEP { class HepRotation; }

namespace mu2e {

  class SimpleConfig;

  void constructExtMonFNALInfrastructure(const VolumeInfo& pixelChillerParent,
                                         const CLHEP::HepRotation& pixelChillerParentRotationInMu2e,
                                         const VolumeInfo& mainParent,
                                         const CLHEP::HepRotation& mainParentRotationInMu2e,
                                         const SimpleConfig& config);

}

#endif /* CONSTRUCTEXTMONFNAL_HH */
