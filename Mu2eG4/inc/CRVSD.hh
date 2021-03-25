#ifndef Mu2eG4_CRVSD_hh
#define Mu2eG4_CRVSD_hh
//
// Define a sensitive detector for CRV Scintillator Bars
//
// Original author KLG
//

// Mu2e includes
#include "Mu2eG4/inc/Mu2eG4SensitiveDetector.hh"

namespace mu2e {

  class CRVSD : public Mu2eG4SensitiveDetector{

  public:

    CRVSD(G4String name, SimpleConfig const & config );
    G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;

  };

} // namespace mu2e

#endif /* Mu2eG4_CRVSD_hh */
