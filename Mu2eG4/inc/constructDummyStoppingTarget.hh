#ifndef constructDummyStoppingTarget_HH
#define constructDummyStoppingTarget_HH
//
// Free function to construct a placeholder for the stopping target.
// Useful for some low detail graphics.
//
// $Id: constructDummyStoppingTarget.hh,v 1.2 2010/12/22 17:37:57 genser Exp $
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

  VolumeInfo constructDummyStoppingTarget( G4LogicalVolume* mother,
                                           double zOff,
                                           const SimpleConfig& config );

}  // end namespace mu2e

#endif
