#ifndef StrawSD_h
#define StrawSD_h 1
//
// Define a sensitive detector for Straws.
// ( Not sure yet if I can use this for both LTracker and TTracker?)
// 
// $Id: StrawSD.hh,v 1.7 2010/09/29 19:37:58 logash Exp $
// $Author: logash $ 
// $Date: 2010/09/29 19:37:58 $
//
// Original author Rob Kutschke
//

// Mu2e includes
#include "Mu2eG4/inc/StepPointG4.hh"
#include "Mu2eG4/inc/EventNumberList.hh"

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
  
  private:

    G4ThreeVector GetTrackerOrigin(const G4TouchableHandle & touchableHandle);

    StepPointG4Collection* _collection;

    // List of events for which to enable debug printout.
    EventNumberList _debugList;
    int _nStrawsPerDevice;
    int _nStrawsPerSector;
    int _TrackerVersion;

    // Limit maximum size of the steps collection
    int _sizeLimit;
    int _currentSize;

  };

} // namespace mu2e

#endif
