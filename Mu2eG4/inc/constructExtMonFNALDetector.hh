// Sam Fine, 2024

#ifndef CONSTRUCTEXTMONFNALDETECTOR_HH
#define CONSTRUCTEXTMONFNALDETECTOR_HH

#include <string>

#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPlaneStack.hh"
#include "Geant4/G4ThreeVector.hh"

namespace CLHEP { class HepRotation; }

namespace mu2e {

  class SimpleConfig;
  class ExtMonFNALMagnet;

  void constructExtMonFNALDetector(const VolumeInfo& parent,
                                   const CLHEP::HepRotation& parentRotationInMu2e,
                                   const SimpleConfig& config
                                  );

  void constructExtMonFNALMagnet(const ExtMonFNALMagnet& mag,
                                 const VolumeInfo& parent,
                                 const std::string& volNameSuffix,
                                 const CLHEP::HepRotation& parentRotationInMu2e,
                                 const SimpleConfig& config
                                 );

  void constructExtMonFNALPlanes(const VolumeInfo& mother,
                                 const ExtMonFNALModule& module,
                                 const ExtMonFNALPlaneStack& stack,
                                 const std::string& volNameSuffix,
                                 const SimpleConfig& config,
                                 bool const forceAuxEdgeVisible,
                                 bool const doSurfaceCheck,
                                 bool const placePV
                                 );

  void constructExtMonFNALModules(const VolumeInfo& mother,
                                  const G4ThreeVector& offset,
                                  unsigned iplane,
                                  const ExtMonFNALModule& module,
                                  const ExtMonFNALPlaneStack& stack,
                                  const std::string& volNameSuffix,
                                  const SimpleConfig& config,
                                  bool const forceAuxEdgeVisible,
                                  bool const doSurfaceCheck,
                                  bool const placePV
                                  );

}


#endif /* CONSTRUCTEXTMONFNALDETECTOR_HH */
