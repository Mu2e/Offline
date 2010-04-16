#ifndef constructDummyTracker_HH
#define constructDummyTracker_HH
//
// Free function to construct a placeholder for the tracker.
// Useful for some low detail graphics.
//
// $Id: constructDummyTracker.hh,v 1.1 2010/04/16 14:45:34 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/16 14:45:34 $
//
// Original author Rob Kutschke
//

#include "Mu2eG4/inc/VolumeInfo.hh"

class G4LogicalVolume;
class G4Material;

namespace mu2e{

  class SimpleConfig;

  VolumeInfo constructDummyTracker( G4LogicalVolume* mother,
                                    double zOff,
                                    const SimpleConfig& config );

}  // end namespace mu2e

#endif
