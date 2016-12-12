#ifndef Mu2eG4_CaloCrystalSD_hh
#define Mu2eG4_CaloCrystalSD_hh
//
// Define a sensitive detector for calorimetric crystals
//
// $Id: CaloCrystalSD.hh,v 1.11 2012/05/29 22:55:26 genser Exp $
// $Author: genser $
// $Date: 2012/05/29 22:55:26 $
//
// Original author Ivan Logashenko
//

#include <map>
#include <vector>

// Mu2e includes
#include "Mu2eG4/inc/EventNumberList.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "Mu2eG4/inc/Mu2eSensitiveDetector.hh"

// Art includes
#include "canvas/Persistency/Provenance/ProductID.h"
#include "art/Framework/Principal/Event.h"

namespace mu2e {

  class CaloCrystalSD : public Mu2eSensitiveDetector{

  public:

    CaloCrystalSD(G4String, const SimpleConfig& config);

    G4bool ProcessHits(G4Step*, G4TouchableHistory*);

  };

} // namespace mu2e

#endif /* Mu2eG4_CaloCrystalSD_hh */
