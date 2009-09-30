#ifndef StrawSD_h
#define StrawSD_h 1
//
// Define a sensitive detector for Straws.
// ( Not sure yet if I can use this for both LTracker and TTracker?)
// 
// $Id: StrawSD.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//

// Mu2e includes
#include "Mu2eG4/inc/StrawG4Hit.hh"

// G4 includes
#include "G4VSensitiveDetector.hh"

class G4Step;
class G4HCofThisEvent;

namespace mu2e {

  class StrawSD : public G4VSensitiveDetector{

  public:
    StrawSD(G4String);
    ~StrawSD();
    
    void Initialize(G4HCofThisEvent*);
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void EndOfEvent(G4HCofThisEvent*);
  
  private:
    StrawG4HitsCollection* _collection;
    
  };

} // namespace mu2e

#endif
