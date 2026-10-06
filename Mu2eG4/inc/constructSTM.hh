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

  void constructSTM(const SimpleConfig& _config);

  // The field of the magnet built by constructSTM().
  // Call from ConstructSDandField(), on every thread.
  void constructSTMMagneticField();

}

#endif /* Mu2eG4_constructSTM_hh */
