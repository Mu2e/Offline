#ifndef Mu2eG4_CaloReadoutSD_hh
#define Mu2eG4_CaloReadoutSD_hh

#include "Mu2eG4/inc/EventNumberList.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "Mu2eG4/inc/Mu2eSensitiveDetector.hh"
#include <map>
#include <vector>


// Art includes
#include "canvas/Persistency/Provenance/ProductID.h"
#include "art/Framework/Principal/Event.h"

namespace mu2e {

  class CaloReadoutSD : public Mu2eSensitiveDetector{

     public:

       CaloReadoutSD(G4String, const SimpleConfig& config);

       G4bool ProcessHits(G4Step*, G4TouchableHistory*);


     private:

       int    _nro;

  };

} 

#endif
