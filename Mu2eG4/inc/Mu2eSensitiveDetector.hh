#ifndef Mu2eG4_Mu2eSensitiveDetector_hh
#define Mu2eG4_Mu2eSensitiveDetector_hh
//
// Defines a generic Mu2e sensitive detector
//
//
// Original author KLG
//

// Mu2e includes
#include "Mu2eG4/inc/EventNumberList.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

// G4 includes
#include "G4VSensitiveDetector.hh"

// Art includes

class G4Step;
class G4HCofThisEvent;

namespace mu2e {

  // Forward declarations in mu2e namespace
  class SimpleConfig;
  class PhysicsProcessInfo;
  class SimParticleHelper;

  class Mu2eSensitiveDetector : public G4VSensitiveDetector{

  public:

    Mu2eSensitiveDetector(G4String const name, SimpleConfig const & config);

    virtual void Initialize(G4HCofThisEvent*) override;

    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;

    virtual void EndOfEvent(G4HCofThisEvent*) override;

    void beforeG4Event(StepPointMCCollection& outputHits,
                       PhysicsProcessInfo & processInfo,
                       const SimParticleHelper& spHelper);

  protected:

    // Non-owning pointer to the  collection into which hits will be added.
    StepPointMCCollection* _collection;

    // Non-ownning pointer and object that returns code describing physics processes.
    PhysicsProcessInfo* _processInfo;

    // Mu2e point of origin
    G4ThreeVector _mu2eOrigin;

    // List of events for which to enable debug printout.
    EventNumberList _debugList;

    // Limit maximum size of the steps collection
    int _sizeLimit;
    int _currentSize;

    // A helper to create pointers to SimParticles
    const SimParticleHelper *_spHelper;
  };

} // namespace mu2e

#endif /* Mu2eG4_Mu2eSensitiveDetector_hh */
