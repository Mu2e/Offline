#ifndef CONSTRUCTEXTMONFNAL_HH
#define CONSTRUCTEXTMONFNAL_HH

namespace CLHEP { class HepRotation; }

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;

  void constructExtMonFNAL(const VolumeInfo& parent,
                           const CLHEP::HepRotation& parentRotationInMu2e,
                           const SimpleConfig& config);

}

#endif /* CONSTRUCTEXTMONFNAL_HH */
