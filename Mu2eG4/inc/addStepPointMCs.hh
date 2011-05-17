#ifndef Mu2eG4_addStepPointMCs_hh
#define Mu2eG4_addStepPointMCs_hh
//
// Add StepPointMC objects to the event.
//
// $Id: addStepPointMCs.hh,v 1.3 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
//
// Original author Rob Kutschke
//

#include "ToyDP/inc/StepPointMCCollection.hh"

class G4Event;

namespace mu2e{

  // Public entry point.
  void addStepPointMCs ( const G4Event *, StepPointMCCollection& hits );

  // Specializations for particular trackers.  Called by the public entry point.
  void addLT( const G4Event* g4event, StepPointMCCollection& hits  );
  void addI ( const G4Event* g4event, StepPointMCCollection& hits  );

}  // end namespace mu2e

#endif /* Mu2eG4_addStepPointMCs_hh */


