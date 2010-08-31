#ifndef AddVirtualDetectorPoints_HH
#define AddVirtualDetectorPoints_HH
//
// Populate output collection for virtual detectors hits
//
// Original author Ivan Logashenko
//

#include "ToyDP/inc/StepPointMCCollection.hh"

class G4Event;

namespace mu2e{

  // Public entry point.
  void addVirtualDetectorPoints ( const G4Event *, StepPointMCCollection& hits );

}  // end namespace mu2e

#endif


