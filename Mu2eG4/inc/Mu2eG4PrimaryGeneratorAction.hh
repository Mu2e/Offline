#ifndef Mu2eG4_Mu2eG4PrimaryGeneratorAction_hh
#define Mu2eG4_Mu2eG4PrimaryGeneratorAction_hh
//
// Give generated tracks to G4 by copying information from a GenParticleCollection.
//
// Original author Rob Kutschke
//

// G4 includes
#include "Geant4/G4VUserPrimaryGeneratorAction.hh"
#include "Geant4/G4ThreeVector.hh"

// Mu2eG4 includes
#include "Mu2eG4/inc/Mu2eG4Config.hh"
#include "DataProducts/inc/PDGCode.hh"

class G4Event;

namespace mu2e {

  class Mu2eG4PerThreadStorage;

  class Mu2eG4PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
  public:

    explicit Mu2eG4PrimaryGeneratorAction(const Mu2eG4Config::Debug& debug,
                                    Mu2eG4PerThreadStorage* pts);

    // This is the interface specified by G4.
    virtual void GeneratePrimaries(G4Event*) override;

  private:

    void addG4Particle(G4Event *event,
                       PDGCode::type pdgId,
                       double excitationEnergy,
                       int floatLevelBaseIndex,
                       const G4ThreeVector& pos,
                       double time,
                       double properTime,
                       const G4ThreeVector& mom);

    int verbosityLevel_;
    // a tuple for ion tests
    std::tuple<int,double,int> testIonToGenerate_;
    // a flag used to enable the above testing
    const bool testPDGIdToGenerate_;

    Mu2eG4PerThreadStorage* perThreadObjects_;

  };

}  // end namespace mu2e
#endif /* Mu2eG4_Mu2eG4PrimaryGeneratorAction_hh */
