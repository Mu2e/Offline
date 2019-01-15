// Original author  K.L.Genser
// A constructor implementing Mu2e specific customizations

// C++ includes

// Framework includes
#include "fhiclcpp/ParameterSet.h"

// Mu2e includes
#include "Mu2eG4/inc/Mu2eG4CustomizationPhysicsConstructor.hh"
#include "Mu2eG4/inc/customizeChargedPionDecay.hh"
#include "Mu2eG4/inc/toggleProcesses.hh"

// G4 includes
#include "G4BuilderType.hh"
#include "G4Threading.hh"
#include "G4EmParameters.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//

namespace mu2e {

  G4_DECLARE_PHYSCONSTR_FACTORY(Mu2eG4CustomizationPhysicsConstructor);

  G4ThreadLocal G4bool Mu2eG4CustomizationPhysicsConstructor::wasActivated = false;

  Mu2eG4CustomizationPhysicsConstructor::
  Mu2eG4CustomizationPhysicsConstructor(const fhicl::ParameterSet* pset)
    : G4VPhysicsConstructor("Mu2eCustomizations",bUnknown),
      pset_(pset)
      // bUnknown disables duplication check for physics builder type
  {}

  // default constructor required by G4PhysicsConstructorFactory; fixme: check, also pset_?
  Mu2eG4CustomizationPhysicsConstructor::
  Mu2eG4CustomizationPhysicsConstructor()
    : G4VPhysicsConstructor("Mu2eCustomizations",bUnknown)
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

    if (pset_->get<int>("debug.diagLevel",0)>0) {
      G4cout << "Mu2eG4CustomizationPhysicsConstructor::"
             << __func__ << " Called" << G4endl;
    }

    if(G4Threading::IsMasterThread()) {

      // G4 does not include pi+ -> e+ nu + cc. Fix that in one of several ways.
      // decay tables are common
      customizeChargedPionDecay(*pset_);

    }

    // Switch off the decay of some particles
    switchDecayOff(*pset_);

    // add conditional and special mu2e processes
    addUserProcesses(*pset_);

    // swap Bertini Cascade with Precompound model in G4MuonMinusCapture if requested
    switchCaptureDModel(*pset_);

    // use more accurate boundary crossing algorithm for muons and hadrons
    // place this late in the setup sequence to avoid a subsequent change
    if (pset_->get<bool>("physics.setMuHadLateralDisplacement",false)) {
      G4EmParameters* params = G4EmParameters::Instance();
      params->SetMuHadLateralDisplacement(true);
    }

  }

}
