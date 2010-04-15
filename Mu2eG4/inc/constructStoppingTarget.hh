#ifndef constructStoppingTarget_HH
#define constructStoppingTarget_HH
//
// Free function to construct the stopping targets.
//
// $Id: constructStoppingTarget.hh,v 1.1 2010/04/15 23:01:39 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/15 23:01:39 $
//
// Original author Peter Shanahan
//
// Notes:

#include "Mu2eG4/inc/VolumeInfo.hh"

class G4LogicalVolume;

namespace mu2e{

  VolumeInfo constructStoppingTarget( G4LogicalVolume* mother,
                                      double zOff );


}  // end namespace mu2e

#endif
