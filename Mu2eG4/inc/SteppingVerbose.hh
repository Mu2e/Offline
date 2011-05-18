#ifndef Mu2eG4_SteppingVerbose_hh
#define Mu2eG4_SteppingVerbose_hh
//
// Verbose version of the stepping action.
//
// $Id: SteppingVerbose.hh,v 1.3 2011/05/18 02:27:17 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:17 $
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
