#ifndef Mu2eG4_TrackerWireSD_hh
#define Mu2eG4_TrackerWireSD_hh
//
// Defines a generic Tracker wire sensitive detector
//
// Original author G. Tassielli
//

//// Mu2e includes
#include "Mu2eG4/inc/Mu2eG4SensitiveDetector.hh"

namespace mu2e {

  class TrackerWireSD : public Mu2eG4SensitiveDetector {

  public:
    TrackerWireSD(G4String, const SimpleConfig& config);
    ~TrackerWireSD();

    virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;

    static void setMu2eDetCenterInWorld(const G4ThreeVector &origin) {
            _mu2eDetCenter = origin;
    }

  protected:

    // Mu2e point of origin
    static G4ThreeVector _mu2eDetCenter;

  };

} // namespace mu2e

#endif /* Mu2eG4_TrackerWireSD_hh */
