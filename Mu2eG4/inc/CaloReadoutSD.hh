#ifndef CaloReadoutSD_h
#define CaloReadoutSD_h 1
//
// Define a sensitive detector for virtual detectors (like G4Beamline)
// 
// Original author Ivan Logashenko
//
#include <map>
#include <vector>

// Mu2e includes
#include "Mu2eG4/inc/StepPointG4.hh"
#include "Mu2eG4/inc/EventNumberList.hh"
#include "Mu2eG4/inc/CaloCrystalSD.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

// G4 includes
#include "G4VSensitiveDetector.hh"

class G4Step;
class G4HCofThisEvent;

namespace mu2e {

  // Forward declarations in mu2e namespace
  class SimpleConfig;

  class CaloReadoutSD : public G4VSensitiveDetector{

  public:
    CaloReadoutSD(G4String, const SimpleConfig& config, CaloCrystalSD *);
    ~CaloReadoutSD();
    
    void Initialize(G4HCofThisEvent*);
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void EndOfEvent(G4HCofThisEvent*);
  
  private:

    CaloCrystalSD *crystalSD;
    int _nro;
    double _minE;
  };

} // namespace mu2e

#endif
