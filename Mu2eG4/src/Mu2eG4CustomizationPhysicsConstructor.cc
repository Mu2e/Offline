// Original author  K.L.Genser
// A constructor implementing Mu2e specific customizations

// C++ includes

// Mu2e includes
#include "Mu2eG4/inc/Mu2eG4CustomizationPhysicsConstructor.hh"
#include "Mu2eG4/inc/customizeChargedPionDecay.hh"
#include "Mu2eG4/inc/toggleProcesses.hh"

// G4 includes
#include "Geant4/G4BuilderType.hh"
#include "Geant4/G4Threading.hh"
#include "Geant4/G4EmParameters.hh"

// factory
#include "Geant4/G4PhysicsConstructorFactory.hh"
//

namespace mu2e {

  G4_DECLARE_PHYSCONSTR_FACTORY(Mu2eG4CustomizationPhysicsConstructor);

  G4ThreadLocal G4bool Mu2eG4CustomizationPhysicsConstructor::wasActivated = false;

  Mu2eG4CustomizationPhysicsConstructor::
  Mu2eG4CustomizationPhysicsConstructor(const Mu2eG4Config::Physics* phys, const Mu2eG4Config::Debug* debug)
    : G4VPhysicsConstructor("Mu2eCustomizations",bUnknown)
    , phys_(phys)
    , debug_(debug)
      // bUnknown disables duplication check for physics builder type
  {}

  // default constructor required by G4PhysicsConstructorFactory; fixme: check, also pset_?
  Mu2eG4CustomizationPhysicsConstructor::
  Mu2eG4CustomizationPhysicsConstructor()
    : G4VPhysicsConstructor("Mu2eCustomizations",bUnknown)
    , phys_(nullptr)
    , debug_(nullptr)
  {}

  Mu2eG4CustomizationPhysicsConstructor::~Mu2eG4CustomizationPhysicsConstructor()
  {}

  void Mu2eG4CustomizationPhysicsConstructor::ConstructParticle()
  {
    // Empty on purpose, for now
  }

  void Mu2eG4CustomizationPhysicsConstructor::ConstructProcess()
  {
    if(wasActivated) { return; }
    wasActivated = true;

    if (debug_->diagLevel()>0) {
      G4cout << "Mu2eG4CustomizationPhysicsConstructor::"
             << __func__ << " Called" << G4endl;
    }

    if(G4Threading::IsMasterThread()) {

      // G4 does not include pi+ -> e+ nu + cc. Fix that in one of several ways.
      // decay tables are common
      customizeChargedPionDecay(*phys_, *debug_);

    }

    // Switch off the decay of some particles
    switchDecayOff(*phys_, *debug_);

    // add conditional and special mu2e processes
    addUserProcesses(*phys_, *debug_);

    // swap Bertini Cascade with Precompound model in G4MuonMinusCapture if requested
    switchCaptureDModel(*phys_, *debug_);

    // use more accurate boundary crossing algorithm for muons and hadrons
    // place this late in the setup sequence to avoid a subsequent change
    if (phys_->setMuHadLateralDisplacement()) {
      G4EmParameters* params = G4EmParameters::Instance();
      params->SetMuHadLateralDisplacement(true);
    }

  }

}
