// Original author  K.L.Genser
// A constructor implementing Mu2e specific customizations
//

#ifndef Mu2eG4CustomizationPhysicsConstructor_h
#define Mu2eG4CustomizationPhysicsConstructor_h 1

#include "G4VPhysicsConstructor.hh"

namespace fhicl { class ParameterSet; }

namespace mu2e {

  class Mu2eG4CustomizationPhysicsConstructor : public G4VPhysicsConstructor
  {
  public:

    explicit Mu2eG4CustomizationPhysicsConstructor(const fhicl::ParameterSet* pset);

    Mu2eG4CustomizationPhysicsConstructor();

    virtual ~Mu2eG4CustomizationPhysicsConstructor();

    virtual void ConstructParticle();

    virtual void ConstructProcess();

  private:

    static G4ThreadLocal G4bool wasActivated;

    // non owning pointer
    // can't be ref due to the default constructor factory requirement
    const fhicl::ParameterSet* pset_;

  };

}
#endif
