//
// Steering routine to G4 call customization routines that must be called
// only after the call to Mu2eG4RunManager::Initialize.
//
// Do not put G4 code in this steering routine.
//
// $Id: postG4InitializeTasks.cc,v 1.4 2012/07/15 22:06:17 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/15 22:06:17 $
//

#include "Mu2eG4/inc/postG4InitializeTasks.hh"

#include "ConfigTools/inc/SimpleConfig.hh"
#include "fhiclcpp/ParameterSet.h"

#include "Mu2eG4/inc/customizeChargedPionDecay.hh"
#include "Mu2eG4/inc/toggleProcesses.hh"
#include "Mu2eG4/inc/setMinimumRangeCut.hh"
#include "Mu2eG4/inc/setParticleCut.hh"

#include "G4VUserPhysicsList.hh"

namespace mu2e{

 template<class Config>
  void postG4InitializeTasks(const Config& config, G4VUserPhysicsList* pL) {

    // G4 does not include pi+ -> e+ nu + cc. Fix that in one of several ways.
    customizeChargedPionDecay(config);

    // Switch off the decay of some particles
    switchDecayOff(config);

    // Set the range cut for e+,e-,gamma (& proton but see below)
    setMinimumRangeCut(config, pL);

    // set the proton production cut which is special
    std::string pName("proton");
    setParticleCut(config, pName, pL);

  }

  template void postG4InitializeTasks(const SimpleConfig&, G4VUserPhysicsList* pL);
  template void postG4InitializeTasks(const fhicl::ParameterSet&, G4VUserPhysicsList* pL);

}  // end namespace mu2e
