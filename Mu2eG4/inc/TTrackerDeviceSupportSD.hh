#ifndef Mu2eG4_TTrackerDeviceSupportSD_hh
#define Mu2eG4_TTrackerDeviceSupportSD_hh
//
// Defines a sensitive detector for TTrackerDeviceSupport
//
// $Id: TTrackerDeviceSupportSD.hh,v 1.2 2012/03/21 15:52:09 gandr Exp $
// $Author: gandr $
// $Date: 2012/03/21 15:52:09 $
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

  class TTrackerDeviceSupportSD : public G4VSensitiveDetector{

  public:
    TTrackerDeviceSupportSD(G4String, const SimpleConfig& config);
    ~TTrackerDeviceSupportSD();

    void Initialize(G4HCofThisEvent*);
    G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*);
    void EndOfEvent(G4HCofThisEvent*);

    void beforeG4Event(StepPointMCCollection& outputHits, 
                       PhysicsProcessInfo & processInfo,
                       art::ProductID const& simID,
                       art::Event const & event );


  private:

    // Non-owning pointer to the  collection into which hits will be added.
    StepPointMCCollection* _collection;

    // Non-ownning pointer and object that returns code describing physics processes.
    PhysicsProcessInfo* _processInfo;

    // Mu2e point of origin
    G4ThreeVector _mu2eOrigin;

    int _TrackerVersion;

    // List of events for which to enable debug printout and verbosityLevel
    EventNumberList _debugList;

    // Limit maximum size of the steps collection
    int _sizeLimit;
    int _currentSize;

    // Information about the SimParticleCollection, needed to instantiate art::Ptr.
    art::ProductID const* _simID;
    art::Event     const* _event;

  };

} // namespace mu2e

#endif /* Mu2eG4_TTrackerDeviceSupportSD_hh */
