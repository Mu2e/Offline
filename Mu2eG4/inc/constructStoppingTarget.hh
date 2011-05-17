#ifndef Mu2eG4_constructStoppingTarget_hh
#define Mu2eG4_constructStoppingTarget_hh
//
// Free function to construct the stopping targets.
//
// $Id: constructStoppingTarget.hh,v 1.4 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
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

#endif /* Mu2eG4_constructStoppingTarget_hh */
