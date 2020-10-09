#ifndef Mu2eG4_checkForOverlaps_hh
#define Mu2eG4_checkForOverlaps_hh
//
// Free function to do Geant4 overlap check
//
//
// Original author KLG
//

class G4VPhysicalVolume;

namespace mu2e{

  class SimpleConfig;

  bool checkForOverlaps( G4VPhysicalVolume* const pv,
                         SimpleConfig const& config, bool verbose=false );

}  // end namespace mu2e

#endif /* Mu2eG4_checkForOverlaps_hh */
