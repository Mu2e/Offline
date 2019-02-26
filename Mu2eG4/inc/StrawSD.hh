#ifndef Mu2eG4_StrawSD_hh
#define Mu2eG4_StrawSD_hh
//
// Define a sensitive detector for Straws.
//
// $Id: StrawSD.hh,v 1.18 2013/12/20 20:09:45 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/12/20 20:09:45 $
//
// Original author Rob Kutschke
//

// Mu2e includes
#include "Mu2eG4/inc/EventNumberList.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "Mu2eG4/inc/Mu2eSensitiveDetector.hh"
#include "TrackerGeom/inc/SupportModel.hh"

// Art includes
#include "canvas/Persistency/Provenance/ProductID.h"
#include "art/Framework/Principal/Event.h"

class G4Step;
class G4HCofThisEvent;

namespace mu2e {

  class StrawSD : public Mu2eSensitiveDetector{

  public:
    StrawSD(const G4String, SimpleConfig const & config);

    G4bool ProcessHits(G4Step*, G4TouchableHistory*);

  private:

    G4ThreeVector GetTrackerOrigin();

    int _nStrawsPerPlane;
    int _nStrawsPerPanel;
    int _TrackerVersion;

    uint16_t _npanels;
    uint16_t _panelsft;
    uint16_t _planesft;

    SupportModel _supportModel;
    int _verbosityLevel;

  };

} // namespace mu2e

#endif /* Mu2eG4_StrawSD_hh */
