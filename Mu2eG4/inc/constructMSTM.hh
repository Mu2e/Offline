#ifndef Mu2eG4_constructMSTM_hh
#define Mu2eG4_constructMSTM_hh
//
// Free function to create the Detector Solenoid
//
//
// Original author KLG
//

namespace mu2e {

  class VolumeInfo;
  class SimpleConfig;

  void constructMSTM(const VolumeInfo& parent,
                     const SimpleConfig& _config
                     );

}

#endif /* Mu2eG4_constructMSTM_hh */
