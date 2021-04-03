#ifndef Mu2eG4_Mu2eG4SteppingVerbose_hh
#define Mu2eG4_Mu2eG4SteppingVerbose_hh
//
// Verbose version of the stepping action.
//
//
// Original author Rob Kutschke
//
// The intial release is just a copy of the G4 Novice N02 example.
//

#include "Geant4/G4SteppingVerbose.hh"

namespace mu2e {
  class Mu2eG4SteppingVerbose : public G4SteppingVerbose {
  public:

    Mu2eG4SteppingVerbose();
    ~Mu2eG4SteppingVerbose();

    void StepInfo();
    void TrackingStarted();

  };

} // end namespace mu2e

#endif /* Mu2eG4_Mu2eG4SteppingVerbose_hh */
