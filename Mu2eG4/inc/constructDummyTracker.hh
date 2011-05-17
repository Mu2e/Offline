#ifndef Mu2eG4_constructDummyTracker_hh
#define Mu2eG4_constructDummyTracker_hh
//
// Free function to construct a placeholder for the tracker.
// Useful for some low detail graphics.
//
// $Id: constructDummyTracker.hh,v 1.3 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
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

#endif /* Mu2eG4_constructDummyTracker_hh */
