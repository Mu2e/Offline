#ifndef Mu2eG4_CaloCrateSD_hh
#define Mu2eG4_CaloCrateSD_hh

// Mu2e includes
#include "Mu2eG4/inc/Mu2eSensitiveDetector.hh"

namespace mu2e {

  class CaloCrateSD : public Mu2eSensitiveDetector{

     public:

       CaloCrateSD(G4String, const SimpleConfig& config);

       G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;

  };

}

#endif
