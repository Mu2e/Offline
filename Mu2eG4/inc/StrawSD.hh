#ifndef Mu2eG4_StrawSD_hh
#define Mu2eG4_StrawSD_hh
//
// Define a sensitive detector for Straws.
//
// Original author Rob Kutschke
//

// Mu2e includes
#include "Mu2eG4/inc/Mu2eG4SensitiveDetector.hh"
#include "TrackerGeom/inc/SupportModel.hh"


namespace mu2e {

  class StrawSD : public Mu2eG4SensitiveDetector{

  public:
    StrawSD(const G4String, SimpleConfig const & config);

    G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;

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
