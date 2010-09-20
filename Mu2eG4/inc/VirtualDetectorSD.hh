#ifndef VirtualDetectorSD_h
#define VirtualDetectorSD_h 1
//
// Define a sensitive detector for virtual detectors (like G4Beamline)
// 
// Original author Ivan Logashenko
//

// Mu2e includes
#include "Mu2eG4/inc/StepPointG4.hh"
#include "Mu2eG4/inc/EventNumberList.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

// G4 includes
#include "G4VSensitiveDetector.hh"

class G4Step;
class G4HCofThisEvent;

namespace mu2e {

  // Forward declarations in mu2e namespace
  class SimpleConfig;

  class VirtualDetectorSD : public G4VSensitiveDetector{

  public:
    VirtualDetectorSD(G4String, const SimpleConfig& config);
    ~VirtualDetectorSD();
    
    void Initialize(G4HCofThisEvent*);
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void EndOfEvent(G4HCofThisEvent*);
  
    static void setMu2eOriginInWorld(const G4ThreeVector &origin) {
      _mu2eOrigin = origin;
    }

    static void setSizeLimit(int sizeLimit) {
      _sizeLimit = sizeLimit;
    }

  private:

    G4ThreeVector GetTrackerOrigin(const G4TouchableHandle & touchableHandle);

    StepPointG4Collection* _collection;

    // List of events for which to enable debug printout.
    EventNumberList _debugList;

    // Mu2e point of origin
    static G4ThreeVector _mu2eOrigin;

    // Limit maximum size of the steps collection
    static int _sizeLimit;
    int _currentSize;
    
  };

} // namespace mu2e

#endif
