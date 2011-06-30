#ifndef Mu2eG4_StrawSD_hh
#define Mu2eG4_StrawSD_hh
//
// Define a sensitive detector for Straws.
//
// $Id: StrawSD.hh,v 1.13 2011/06/30 04:55:13 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/30 04:55:13 $
//
// Original author Rob Kutschke
//

// Mu2e includes
#include "Mu2eG4/inc/EventNumberList.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

// G4 includes
#include "G4VSensitiveDetector.hh"

// Art includes
#include "art/Persistency/Provenance/ProductID.h"
#include "art/Persistency/Common/EDProductGetter.h"

class G4Step;
class G4HCofThisEvent;

namespace mu2e {

  // Forward declarations in mu2e namespace
  class SimpleConfig;
  class PhysicsProcessInfo;

  class StrawSD : public G4VSensitiveDetector{

  public:
    StrawSD(G4String, const SimpleConfig& config);
    ~StrawSD();

    void Initialize(G4HCofThisEvent*);
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void EndOfEvent(G4HCofThisEvent*);

    void beforeG4Event(StepPointMCCollection& outputHits, 
                       PhysicsProcessInfo & processInfo,
                       art::ProductID const& simID, 
                       art::EDProductGetter const* productGetter
                       );

  private:

    G4ThreeVector GetTrackerOrigin(const G4TouchableHandle & touchableHandle);

    // Non-owning pointer to the  collection into which hits will be added.
    StepPointMCCollection* _collection;

    // Non-ownning pointer and object that returns code describing physics processes.
    PhysicsProcessInfo* _processInfo;

    int _nStrawsPerDevice;
    int _nStrawsPerSector;
    int _TrackerVersion;

    // List of events for which to enable debug printout.
    EventNumberList _debugList;

    // Limit maximum size of the steps collection
    int _sizeLimit;
    int _currentSize;

    // Information about the SimParticleCollection, needed to instantiate art::Ptr.
    art::ProductID const *      _simID;
    art::EDProductGetter const* _productGetter;


  };

} // namespace mu2e

#endif /* Mu2eG4_StrawSD_hh */
