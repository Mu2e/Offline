#ifndef Mu2eG4_TrackerPlaneSupportSD_hh
#define Mu2eG4_TrackerPlaneSupportSD_hh
//
// Defines a sensitive detector for TrackerPlaneSupport
//
//
// Original author KLG
//

// Mu2e includes
#include "Offline/Mu2eG4/inc/EventNumberList.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/Mu2eG4/inc/Mu2eG4SensitiveDetector.hh"

// Art includes
#include "canvas/Persistency/Provenance/ProductID.h"
#include "art/Framework/Principal/Event.h"

namespace mu2e {

  class TrackerPlaneSupportSD : public Mu2eG4SensitiveDetector{

  public:

    TrackerPlaneSupportSD(G4String, const SimpleConfig& config);

    G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*);


  private:

    int _TrackerVersion;

  };

} // namespace mu2e

#endif /* Mu2eG4_TrackerPlaneSupportSD_hh */
