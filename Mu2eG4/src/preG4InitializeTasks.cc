//
// Steering routine to G4 call customization routines that must be called
// only after the call to Mu2eG4RunManager::Initialize.
//
// Do not put G4 code in this steering routine.
//


#include "Mu2eG4/inc/preG4InitializeTasks.hh"
//#include "Mu2eG4/inc/setBirksConstant.hh"

namespace mu2e{

  void preG4InitializeTasks(const Mu2eG4Config::Top& config) {

    // Modify Birks Constant for some materials
    // doing it here does not work for Mu2e materials; moving it to ConstructAterials
    //    setBirksConstant(config);

  }

}  // end namespace mu2e
