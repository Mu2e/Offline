#ifndef constructDummyStoppingTarget_HH
#define constructDummyStoppingTarget_HH
//
// Free function to construct a placeholder for the stopping target.
// Useful for some low detail graphics.
//
// $Id: constructDummyStoppingTarget.hh,v 1.1 2010/04/16 14:45:34 kutschke Exp $
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

  VolumeInfo constructDummyStoppingTarget( G4LogicalVolume* mother,
                                           double zOff,
                                           const SimpleConfig& config );

}  // end namespace mu2e

#endif
