// Andrei Gaponenko, 2012

#ifndef CONSTRUCTEXTMONFNAL_HH
#define CONSTRUCTEXTMONFNAL_HH

#include <string>

#include "Mu2eG4Helper/inc/VolumeInfo.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPlaneStack.hh"
#include "Geant4/G4ThreeVector.hh"

namespace CLHEP { class HepRotation; }

namespace mu2e {

  class SimpleConfig;
  class ExtMonFNALMagnet;
  class ExtMonFNALMuonID;

  void constructExtMonFNAL(const VolumeInfo& collimator1Parent,
                           const CLHEP::HepRotation& collimator1ParentRotationInMu2e,
                           const VolumeInfo& mainParent,
                           const CLHEP::HepRotation& mainParentRotationInMu2e,
                           const SimpleConfig& config);

  void constructExtMonFNALBuilding(const VolumeInfo& collimator1Parent,
                                   const CLHEP::HepRotation& collimator1ParentRotationInMu2e,
                                   const VolumeInfo& mainParent,
                                   const CLHEP::HepRotation& mainParentRotationInMu2e,
                                   const SimpleConfig& config);

  void constructExtMonFNALDetector(const VolumeInfo& room, const SimpleConfig& config);

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

  void constructExtMonFNALMuonID(const ExtMonFNALModule& module,
                                 const ExtMonFNALMuonID& muid,
                                 const std::string& volNameSuffix,
                                 const VolumeInfo& parent,
                                 const CLHEP::HepRotation& parentRotationInMu2e,
                                 const SimpleConfig& config
                                 );

}

#endif /* CONSTRUCTEXTMONFNAL_HH */
