#ifndef Mu2eG4BiasedRDPhysics_h
#define Mu2eG4BiasedRDPhysics_h 1

#include "Geant4/G4VPhysicsConstructor.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4Config.hh"

namespace mu2e{

  class Mu2eG4BiasedRDPhysics : public G4VPhysicsConstructor
  {
    public:
      Mu2eG4BiasedRDPhysics();
      Mu2eG4BiasedRDPhysics(const Mu2eG4Config::Physics* phys, G4int verbose=0);

     ~Mu2eG4BiasedRDPhysics() override = default;

      void ConstructParticle() override;
      void ConstructProcess()  override;

   private:

    // non owning pointer
    // can't be ref due to the default constructor factory requirement
    const Mu2eG4Config::Physics* phys_;
    G4int                        verbose_;

  };
}


#endif
