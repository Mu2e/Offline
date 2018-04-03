#ifndef Mu2eG4_constructSTM_hh
#define Mu2eG4_constructSTM_hh
//
// Free function to create Stopping Target Monitor
//
//
// Author: Anthony Palladino
//

namespace mu2e {

  class SimpleConfig;
  class SensitiveDetectorHelper;

  void constructSTM(const SimpleConfig& _config,
                    const SensitiveDetectorHelper& sdHelper
                    );

}

#endif /* Mu2eG4_constructSTM_hh */
