#ifndef SteppingVerbose_h
#define SteppingVerbose_h 1
//
// Verbose version of the stepping action.
// 
// $Id: SteppingVerbose.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
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

#endif
