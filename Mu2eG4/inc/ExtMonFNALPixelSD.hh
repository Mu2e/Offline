#ifndef Mu2eG4_ExtMonFNALPixelSD_hh
#define Mu2eG4_ExtMonFNALPixelSD_hh
//
// $Id: ExtMonFNALPixelSD.hh,v 1.1 2012/08/23 23:36:14 gandr Exp $
// $Author: gandr $
// $Date: 2012/08/23 23:36:14 $
//
// Define a sensitive detector for the FNAL extinction monitor
// Andrei Gaponenko, 2012
//

#include "G4VSensitiveDetector.hh"

#include "art/Persistency/Provenance/ProductID.h"
#include "art/Framework/Principal/Event.h"

namespace mu2e {

  class SimpleConfig;
  class ExtMonFNALSimHitCollection;

  // This class should not derive from Mu2eSensitiveDetector because the latter
  // presumes a wrong type for the hit collection.
  class ExtMonFNALPixelSD : public G4VSensitiveDetector {

    // Non-owning pointer to the  collection into which hits will be added.
    ExtMonFNALSimHitCollection* _collection;

    // Limit maximum size of the steps collection
    unsigned _sizeLimit;

    // Information about the SimParticleCollection, needed to instantiate art::Ptr.
    art::ProductID const * _simID;
    art::Event     const * _event;

  public:

    explicit ExtMonFNALPixelSD(const SimpleConfig& config);

    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

    virtual void EndOfEvent(G4HCofThisEvent*);

    void beforeG4Event(ExtMonFNALSimHitCollection *outputHits,
                       const art::ProductID& simID,
                       const art::Event& event);
  };

} // namespace mu2e

#endif /* Mu2eG4_ExtMonFNALPixelSD_hh */
