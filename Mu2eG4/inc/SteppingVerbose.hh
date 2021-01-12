#ifndef Mu2eG4_SteppingVerbose_hh
#define Mu2eG4_SteppingVerbose_hh
//
// Verbose version of the stepping action.
//
//
// Original author Rob Kutschke
//
// The intial release is just a copy of the G4 Novice N02 example.
//

#include "G4SteppingVerbose.hh"

namespace mu2e {
  class SteppingVerbose : public G4SteppingVerbose {
  public:

    SteppingVerbose();
    ~SteppingVerbose();

    void StepInfo();
    void TrackingStarted();

  };

} // end namespace mu2e

#endif /* Mu2eG4_SteppingVerbose_hh */
