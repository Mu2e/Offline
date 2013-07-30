#ifndef Mu2eG4_ExtMonFNALPixelSD_hh
#define Mu2eG4_ExtMonFNALPixelSD_hh
//
// $Id: ExtMonFNALPixelSD.hh,v 1.2 2013/07/30 18:45:00 wieschie Exp $
// $Author: wieschie $
// $Date: 2013/07/30 18:45:00 $
//
// Define a sensitive detector for the FNAL extinction monitor
// Andrei Gaponenko, 2012
//

#include "G4VSensitiveDetector.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
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
    
    const mu2e::ExtMonFNAL::ExtMon *_extmon;
    
  public:

    ExtMonFNALPixelSD(const SimpleConfig& config, const mu2e::ExtMonFNAL::ExtMon& extmon);

    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

    virtual void EndOfEvent(G4HCofThisEvent*);

    void beforeG4Event(ExtMonFNALSimHitCollection *outputHits,
                       const art::ProductID& simID,
                       const art::Event& event);
  };

} // namespace mu2e

#endif /* Mu2eG4_ExtMonFNALPixelSD_hh */
