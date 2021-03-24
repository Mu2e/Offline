#ifndef Mu2eG4_CaloReadoutCardSD_hh
#define Mu2eG4_CaloReadoutCardSD_hh

// Mu2e includes
#include "Mu2eG4/inc/Mu2eG4SensitiveDetector.hh"

namespace mu2e {

  class CaloReadoutCardSD : public Mu2eG4SensitiveDetector{

     public:

       CaloReadoutCardSD(G4String, const SimpleConfig& config);
       G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;

     private:

       int  _nro;

  };

}

#endif
