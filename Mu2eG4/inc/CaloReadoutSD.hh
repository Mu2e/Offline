#ifndef Mu2eG4_CaloReadoutSD_hh
#define Mu2eG4_CaloReadoutSD_hh

#include "Mu2eG4/inc/Mu2eG4SensitiveDetector.hh"


namespace mu2e {

  class CaloReadoutSD : public Mu2eG4SensitiveDetector{

     public:

       CaloReadoutSD(G4String, const SimpleConfig& config);

       G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;


     private:

       int    _nro;

  };

}

#endif
