//
// Steering routine to G4 call customization routines that must be called
// only after the call to Mu2eG4RunManager::Initialize.
//
// Do not put G4 code in this steering routine.
//
// $Id: postG4InitializationTasks.cc,v 1.1 2012/06/04 19:28:01 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/06/04 19:28:01 $
//

#include "Mu2eG4/inc/postG4InitializeTasks.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

#include "Mu2eG4/inc/toggleProcesses.hh"
#include "Mu2eG4/inc/setMinimumRangeCut.hh"

namespace mu2e{

  void postG4InitializeTasks( SimpleConfig const& config ){

    // Switch off the decay of some particles
    switchDecayOff(config);

    // Add user processes
    addUserProcesses(config);

    // If requested, change the minimum range cut.
    setMinimumRangeCut(config);

  }

}  // end namespace mu2e
