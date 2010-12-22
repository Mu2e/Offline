#ifndef constructStoppingTarget_HH
#define constructStoppingTarget_HH
//
// Free function to construct the stopping targets.
//
// $Id: constructStoppingTarget.hh,v 1.3 2010/12/22 17:37:57 genser Exp $
// $Author: genser $
// $Date: 2010/12/22 17:37:57 $
//
// Original author Peter Shanahan
//
// Notes:

#include "G4Helper/inc/VolumeInfo.hh"

class G4LogicalVolume;

namespace mu2e{

  class SimpleConfig;

  VolumeInfo constructStoppingTarget( G4LogicalVolume* mother,
                                      double zOff,
                                      SimpleConfig const& config );


}  // end namespace mu2e

#endif
