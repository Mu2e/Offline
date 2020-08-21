#ifndef Mu2eG4_constructCRV_hh
#define Mu2eG4_constructCRV_hh
//
// Free function to create CRV aka Scintillator Shield in CosmicRayShield
//
//
// Original author KLG
//

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;

  void constructCRV(VolumeInfo   const & parent,
                    SimpleConfig const & _config
                    );

}

#endif /* Mu2eG4_constructCRV_hh */
