#ifndef Mu2eG4_Mu2eSensitiveDetector_hh
#define Mu2eG4_Mu2eSensitiveDetector_hh
//
// Defines a generic Mu2e sensitive detector
//
// $Id: Mu2eSensitiveDetector.hh,v 1.2 2013/08/28 05:58:17 gandr Exp $
// $Author: gandr $
// $Date: 2013/08/28 05:58:17 $
//
// Original author KLG
//

// Mu2e includes
#include "Mu2eG4/inc/EventNumberList.hh"
#include "Mu2eG4/inc/SimParticleHelper.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

// G4 includes
#include "G4VSensitiveDetector.hh"

// Art includes
#include "canvas/Persistency/Provenance/ProductID.h"
#include "art/Framework/Principal/Event.h"

class G4Step;
class G4HCofThisEvent;

namespace mu2e {

  // Forward declarations in mu2e namespace
  class SimpleConfig;
  class PhysicsProcessInfo;

  class Mu2eSensitiveDetector : public G4VSensitiveDetector{

  public:

    Mu2eSensitiveDetector(G4String const name, SimpleConfig const & config);

    virtual void Initialize(G4HCofThisEvent*);

    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

    virtual void EndOfEvent(G4HCofThisEvent*);

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
