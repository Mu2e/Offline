#ifndef Mu2eG4_TrackerPlaneSupportSD_hh
#define Mu2eG4_TrackerPlaneSupportSD_hh
//
// Defines a sensitive detector for TrackerPlaneSupport
//
// $Id: TrackerPlaneSupportSD.hh,v 1.3 2012/05/29 22:57:44 genser Exp $
// $Author: genser $
// $Date: 2012/05/29 22:57:44 $
//
// Original author KLG
//

// Mu2e includes
#include "Mu2eG4/inc/EventNumberList.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "Mu2eG4/inc/Mu2eSensitiveDetector.hh"

// Art includes
#include "canvas/Persistency/Provenance/ProductID.h"
#include "art/Framework/Principal/Event.h"

namespace mu2e {

  class TrackerPlaneSupportSD : public Mu2eSensitiveDetector{

  public:

    TrackerPlaneSupportSD(G4String, const SimpleConfig& config);

    G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*);


  private:

    int _TrackerVersion;

  };

} // namespace mu2e

#endif /* Mu2eG4_TrackerPlaneSupportSD_hh */
