#ifndef Mu2eG4_ExtMonFNALPixelSD_hh
#define Mu2eG4_ExtMonFNALPixelSD_hh
//
//
// Define a sensitive detector for the FNAL extinction monitor
// Andrei Gaponenko, 2012
//

#include "G4VSensitiveDetector.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "canvas/Persistency/Provenance/ProductID.h"
#include "art/Framework/Principal/Event.h"

#include "MCDataProducts/inc/ExtMonFNALSimHitCollection.hh"

namespace mu2e {

  class SimpleConfig;
  class SimParticleHelper;

  // This class should not derive from Mu2eSensitiveDetector because the latter
  // presumes a wrong type for the hit collection.
  class ExtMonFNALPixelSD : public G4VSensitiveDetector {

    // Non-owning pointer to the  collection into which hits will be added.
    ExtMonFNALSimHitCollection* _collection;

    // Limit maximum size of the steps collection
    unsigned _sizeLimit;

    // Information about the SimParticleCollection, needed to instantiate art::Ptr.
    const SimParticleHelper *_spHelper;

    const mu2e::ExtMonFNAL::ExtMon *_extmon;

  public:

    ExtMonFNALPixelSD(const SimpleConfig& config, const mu2e::ExtMonFNAL::ExtMon& extmon);

    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

    virtual void EndOfEvent(G4HCofThisEvent*);

    void beforeG4Event(ExtMonFNALSimHitCollection *outputHits,
                       const SimParticleHelper& spHelper);
  };

} // namespace mu2e

#endif /* Mu2eG4_ExtMonFNALPixelSD_hh */
