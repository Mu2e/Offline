//
// Steering routine to G4 call customization routines that must be called
// only after the call to Mu2eG4RunManager::Initialize.
//
// Do not put G4 code in this steering routine.
//
// $Id: preG4InitializeTasks.cc,v 1.4 2012/11/04 22:06:17 genser Exp $
// $Author: genser $
// $Date: 2015/11/04 22:06:17 $
//

#include "Mu2eG4/inc/preG4InitializeTasks.hh"

#include "fhiclcpp/ParameterSet.h"

#include "Mu2eG4/inc/setBirksConstant.hh"

namespace mu2e{

  void preG4InitializeTasks(const fhicl::ParameterSet& pset) {

    // Modify Birks Constant for some materials
    setBirksConstant(pset);

  }

}  // end namespace mu2e
