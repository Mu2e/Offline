#ifndef AddStepPointMCs_HH
#define AddStepPointMCs_HH
//
// Add StepPointMC objects to the event.
//
// $Id: addStepPointMCs.hh,v 1.1 2010/04/07 23:19:57 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/07 23:19:57 $
//
// Original author Rob Kutschke
//

#include "ToyDP/inc/StepPointMCCollection.hh"

class G4Event;

namespace mu2e{

  // Public entry point.
  void addStepPointMCs ( const G4Event *, StepPointMCCollection& hits );

  // Specializations for particular trackers.  Called by the public entry point.
  void addL( const G4Event* g4event, StepPointMCCollection& hits  );
  void addI( const G4Event* g4event, StepPointMCCollection& hits  );

}  // end namespace mu2e

#endif


