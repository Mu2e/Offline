#ifndef Mu2eG4_checkForOverlaps_hh
#define Mu2eG4_checkForOverlaps_hh
//
// Free function to do Geant4 overlap check
//
// $Id: checkForOverlaps.hh,v 1.1 2013/08/30 15:52:41 genser Exp $
// $Author: genser $
// $Date: 2013/08/30 15:52:41 $
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
