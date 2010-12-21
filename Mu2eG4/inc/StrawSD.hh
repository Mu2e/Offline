#ifndef StrawSD_h
#define StrawSD_h 1
//
// Define a sensitive detector for Straws.
// 
// $Id: StrawSD.hh,v 1.8 2010/12/21 21:47:51 genser Exp $
// $Author: genser $ 
// $Date: 2010/12/21 21:47:51 $
//
// Original author Rob Kutschke
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

  class StrawSD : public G4VSensitiveDetector{

  public:
    StrawSD(G4String, const SimpleConfig& config);
    ~StrawSD();
    
    void Initialize(G4HCofThisEvent*);
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void EndOfEvent(G4HCofThisEvent*);
  
    void beforeG4Event(StepPointMCCollection& outputHits);

  private:

    G4ThreeVector GetTrackerOrigin(const G4TouchableHandle & touchableHandle);

    StepPointMCCollection* _collection;

    int _nStrawsPerDevice;
    int _nStrawsPerSector;
    int _TrackerVersion;

    // List of events for which to enable debug printout.
    EventNumberList _debugList;

    // Limit maximum size of the steps collection
    int _sizeLimit;
    int _currentSize;

  };

} // namespace mu2e

#endif
