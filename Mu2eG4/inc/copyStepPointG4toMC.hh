#ifndef Mu2eG4_copyStepPointG4toMC_hh
#define Mu2eG4_copyStepPointG4toMC_hh
//
// Populate output StepPointMC collection from StepPointG4
//
// Original author Ivan Logashenko
//
#include <string>

#include "ToyDP/inc/StepPointMCCollection.hh"

class G4Event;

namespace mu2e{

  // Public entry point.
  void copyStepPointG4toMC ( const G4Event *, const std::string, StepPointMCCollection& );

}  // end namespace mu2e

#endif /* Mu2eG4_copyStepPointG4toMC_hh */


