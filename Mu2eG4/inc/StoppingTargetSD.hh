#ifndef Mu2eG4_StoppingTargetSD_hh
#define Mu2eG4_StoppingTargetSD_hh
//
// Define a sensitive detector for stopping target
//
// Original author Ivan Logashenko
//

// Mu2e includes
#include "Mu2eG4/inc/EventNumberList.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"

// G4 includes
#include "G4VSensitiveDetector.hh"

class G4Step;
class G4HCofThisEvent;

namespace mu2e {

  // Forward declarations in mu2e namespace
  class SimpleConfig;

  class StoppingTargetSD : public G4VSensitiveDetector{

  public:
    StoppingTargetSD(G4String, const SimpleConfig& config);
    ~StoppingTargetSD();

    void Initialize(G4HCofThisEvent*);
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void EndOfEvent(G4HCofThisEvent*);

    void beforeG4Event(StepPointMCCollection& outputHits);

    static void setMu2eOriginInWorld(const G4ThreeVector &origin) {
      _mu2eOrigin = origin;
    }

  private:

    StepPointMCCollection* _collection;

    // Mu2e point of origin
    static G4ThreeVector _mu2eOrigin;

    // List of events for which to enable debug printout.
    EventNumberList _debugList;

    // Limit maximum size of the steps collection
    int _sizeLimit;
    int _currentSize;

  };

} // namespace mu2e

#endif /* Mu2eG4_StoppingTargetSD_hh */
