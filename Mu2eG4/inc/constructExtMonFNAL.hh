// Andrei Gaponenko, 2012

#ifndef CONSTRUCTEXTMONFNAL_HH
#define CONSTRUCTEXTMONFNAL_HH

#include <string>

#include "G4Helper/inc/VolumeInfo.hh"

namespace CLHEP { class HepRotation; }

namespace mu2e {

  class SimpleConfig;
  class ExtMonFNALMagnet;

  void constructExtMonFNAL(const VolumeInfo& collimator1Parent,
                           const CLHEP::HepRotation& collimator1ParentRotationInMu2e,
                           const VolumeInfo& mainParent,
                           const CLHEP::HepRotation& mainParentRotationInMu2e,
                           const SimpleConfig& config);

  VolumeInfo constructExtMonFNALBuilding(const VolumeInfo& collimator1Parent,
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

}

#endif /* CONSTRUCTEXTMONFNAL_HH */
