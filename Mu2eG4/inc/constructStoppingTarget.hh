#ifndef constructStoppingTarget_HH
#define constructStoppingTarget_HH
//
// Free function to construct the stopping targets.
//
// $Id: constructStoppingTarget.hh,v 1.2 2010/08/31 16:13:15 genser Exp $
// $Author: genser $
// $Date: 2010/08/31 16:13:15 $
//
// Original author Peter Shanahan
//
// Notes:

#include "Mu2eG4/inc/VolumeInfo.hh"

class G4LogicalVolume;

namespace mu2e{

  class SimpleConfig;

  VolumeInfo constructStoppingTarget( G4LogicalVolume* mother,
                                      double zOff,
                                      SimpleConfig const& config );


}  // end namespace mu2e

#endif
