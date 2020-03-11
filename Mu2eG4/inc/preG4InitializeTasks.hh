#ifndef Mu2eG4_preG4InitializeTasks_hh
#define Mu2eG4_preG4InitializeTasks_hh
//
// Steering routine to G4 call customization routines that must be called
// prior to the call to Mu2eG4RunManager::Initialize.
//
// $Id: preG4InitializeTasks.hh,v 1.1 2012/06/04 19:28:01 genser Exp $
// $Author: genser $
// $Date: 2012/06/04 19:28:01 $
//

#include "Mu2eG4/inc/Mu2eG4Config.hh"

namespace mu2e{

  void preG4InitializeTasks(const Mu2eG4Config::Physics& phys, const Mu2eG4Config::Debug& debug);

}  // end namespace mu2e

#endif /* Mu2eG4_preG4InitializeTasks_hh */
