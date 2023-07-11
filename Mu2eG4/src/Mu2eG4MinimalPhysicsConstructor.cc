//
// A PhysicsConstructor that implements the Minimal physics
//
// Original author KLG inspired by Rob Kutschke's Minimal Physics list
//
// Mu2e includes
#include "Offline/Mu2eG4/inc/Mu2eG4MinimalPhysicsConstructor.hh"
// G4 includes
#include "Geant4/G4ParticleTypes.hh"
#include "Geant4/G4PhysicsConstructorFactory.hh"

namespace mu2e {

  G4_DECLARE_PHYSCONSTR_FACTORY(Mu2eG4MinimalPhysicsConstructor);

  // The unused second parameter (physics type) defaults to bUnknown =0
  Mu2eG4MinimalPhysicsConstructor::Mu2eG4MinimalPhysicsConstructor():
    G4VPhysicsConstructor("Mu2eG4MinimalPhysicsConstructor"){
  }

  void Mu2eG4MinimalPhysicsConstructor::ConstructParticle(){

    G4ChargedGeantino::Definition();
    G4Electron::Definition();
    G4Positron::Definition();
    G4MuonPlus::Definition();
    G4MuonMinus::Definition();
    G4Gamma::Definition();
    G4Proton::Definition();
    G4AntiProton::Definition();
    G4GenericIon::Definition();

   }

  void Mu2eG4MinimalPhysicsConstructor::ConstructProcess(){

    // planing to add the use Mu2eG4Config in future to be able to
    // control the printout
    // if (debug_->diagLevel()>0) {
    //   G4cout << "Mu2eG4MinimalPhysicsConstructor::"
    //          << __func__ << " Called" << G4endl;
    // }

  }


}  // end namespace mu2e
