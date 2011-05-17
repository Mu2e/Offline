#ifndef Mu2eG4_SteppingVerbose_hh
#define Mu2eG4_SteppingVerbose_hh
//
// Verbose version of the stepping action.
// 
// $Id: SteppingVerbose.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $ 
// $Date: 2011/05/17 15:41:36 $
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
