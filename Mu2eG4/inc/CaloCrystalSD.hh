#ifndef Mu2eG4_CaloCrystalSD_hh
#define Mu2eG4_CaloCrystalSD_hh
//
// Define a sensitive detector for calorimetric crystals
//
// $Id: CaloCrystalSD.hh,v 1.5 2011/05/18 02:27:17 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:17 $
//
// Original author Ivan Logashenko
//

#include <map>
#include <vector>

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

  class CaloCrystalSD : public G4VSensitiveDetector{

  public:
    CaloCrystalSD(G4String, const SimpleConfig& config);
    ~CaloCrystalSD();

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

#endif /* Mu2eG4_CaloCrystalSD_hh */
