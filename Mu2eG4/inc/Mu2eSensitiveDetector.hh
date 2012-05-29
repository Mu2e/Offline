#ifndef Mu2eG4_Mu2eSensitiveDetector_hh
#define Mu2eG4_Mu2eSensitiveDetector_hh
//
// Defines a generic Mu2e sensitive detector
//
// $Id: Mu2eSensitiveDetector.hh,v 1.1 2012/05/29 22:52:49 genser Exp $
// $Author: genser $
// $Date: 2012/05/29 22:52:49 $
//
// Original author KLG
//

// Mu2e includes
#include "Mu2eG4/inc/EventNumberList.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

// G4 includes
#include "G4VSensitiveDetector.hh"

// Art includes
#include "art/Persistency/Provenance/ProductID.h"
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
                       art::ProductID const& simID,
                       art::Event const & event );

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

    // Information about the SimParticleCollection, needed to instantiate art::Ptr.
    art::ProductID const * _simID;
    art::Event     const * _event;

  };

} // namespace mu2e

#endif /* Mu2eG4_Mu2eSensitiveDetector_hh */
