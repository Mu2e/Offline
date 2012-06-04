//
// Steering routine to G4 call customization routines that must be called
// only after the call to Mu2eG4RunManager::Initialize.
//
// Do not put G4 code in this steering routine.
//
// $Id: postG4InitializeTasks.cc,v 1.2 2012/06/04 23:50:50 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/06/04 23:50:50 $
//

#include "Mu2eG4/inc/postG4InitializeTasks.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

#include "Mu2eG4/inc/toggleProcesses.hh"
#include "Mu2eG4/inc/setMinimumRangeCut.hh"

#include "Mu2eG4/inc/checkMSCmodel.hh"

namespace mu2e{

  void postG4InitializeTasks( SimpleConfig const& config ){

    // Switch off the decay of some particles
    switchDecayOff(config);

    // Add user processes
    addUserProcesses(config);

    // If requested, change the minimum range cut.
    setMinimumRangeCut(config);

    // If the ITracker is used, check the geant4 Multiple Scattering Model selected.
    checkMSCmodel(config);

  }

}  // end namespace mu2e
