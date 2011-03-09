#ifndef CRSScintillatorBarSD_h
#define CRSScintillatorBarSD_h
//
// Define a sensitive detector for 
// 
// $Id: CRSScintillatorBarSD.hh,v 1.1 2011/03/09 19:25:27 genser Exp $
// $Author: genser $ 
// $Date: 2011/03/09 19:25:27 $
// 
// Original author KLG
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

  class CRSScintillatorBarSD : public G4VSensitiveDetector{

  public:
    CRSScintillatorBarSD(G4String, const SimpleConfig& config);
    ~CRSScintillatorBarSD();
    
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

#endif
