#ifndef Mu2eG4_TrackerWireSD_hh
#define Mu2eG4_TrackerWireSD_hh
//
// Defines a generic Tracker wire sensitive detector
//
// $Id: TrackerWireSD.hh,v 1.2 2012/12/04 00:51:28 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/12/04 00:51:28 $
//
// Original author G. Tassielli
//

//// Mu2e includes
#include "Mu2eG4/inc/Mu2eSensitiveDetector.hh"

namespace mu2e {

  class TrackerWireSD : public Mu2eSensitiveDetector {

  public:
    TrackerWireSD(G4String, const SimpleConfig& config);
    ~TrackerWireSD();

    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

    static void setMu2eDetCenterInWorld(const G4ThreeVector &origin) {
            _mu2eDetCenter = origin;
    }

  protected:

    // Mu2e point of origin
    static G4ThreeVector _mu2eDetCenter;

  };

} // namespace mu2e

#endif /* Mu2eG4_TrackerWireSD_hh */
