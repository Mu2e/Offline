#ifndef Mu2eG4_CaloCrystalSD_hh
#define Mu2eG4_CaloCrystalSD_hh
//
// Define a sensitive detector for calorimetric crystals
//
// Original author Ivan Logashenko
//

// Mu2e includes
#include "Mu2eG4/inc/Mu2eG4SensitiveDetector.hh"

namespace mu2e {

  class CaloCrystalSD : public Mu2eG4SensitiveDetector{

  public:

    CaloCrystalSD(G4String, const SimpleConfig& config);

    G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;

  };

} // namespace mu2e

#endif /* Mu2eG4_CaloCrystalSD_hh */
