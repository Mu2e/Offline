#ifndef constructDummyTracker_HH
#define constructDummyTracker_HH
//
// Free function to construct a placeholder for the tracker.
// Useful for some low detail graphics.
//
// $Id: constructDummyTracker.hh,v 1.2 2010/12/22 17:37:57 genser Exp $
// $Author: genser $
// $Date: 2010/12/22 17:37:57 $
//
// Original author Rob Kutschke
//

#include "G4Helper/inc/VolumeInfo.hh"

class G4LogicalVolume;
class G4Material;

namespace mu2e{

  class SimpleConfig;

  VolumeInfo constructDummyTracker( G4LogicalVolume* mother,
                                    double zOff,
                                    const SimpleConfig& config );

}  // end namespace mu2e

#endif
