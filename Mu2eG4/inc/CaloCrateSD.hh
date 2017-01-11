#ifndef Mu2eG4_CaloCrateSD_hh
#define Mu2eG4_CaloCrateSD_hh

// Mu2e includes
#include "Mu2eG4/inc/EventNumberList.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "Mu2eG4/inc/Mu2eSensitiveDetector.hh"

// Art includes
#include "canvas/Persistency/Provenance/ProductID.h"
#include "art/Framework/Principal/Event.h"

namespace mu2e {

  class CaloCrateSD : public Mu2eSensitiveDetector{

     public:

       CaloCrateSD(G4String, const SimpleConfig& config);
       G4bool ProcessHits(G4Step*, G4TouchableHistory*);

  };

}

#endif
