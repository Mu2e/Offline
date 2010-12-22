#ifndef constructTTracker_HH
#define constructTTracker_HH
//
// Free functions to construct various versions of the TTracker.
//
// $Id: constructTTracker.hh,v 1.4 2010/12/22 17:37:57 genser Exp $
// $Author: genser $
// $Date: 2010/12/22 17:37:57 $
//
// Original author Rob Kutschke
//

#include "G4Helper/inc/VolumeInfo.hh"

class G4LogicalVolume;

namespace mu2e{

  class SimpleConfig;

  VolumeInfo constructTTrackerv1( G4LogicalVolume* mother,
                                  double zOff,
                                  SimpleConfig const& config );

  VolumeInfo constructTTrackerv2( G4LogicalVolume* mother,
                                  double zOff,
                                  SimpleConfig const& config );

  VolumeInfo constructTTrackerv3( G4LogicalVolume* mother,
                                  double zOff,
                                  SimpleConfig const& config );

}  // end namespace mu2e

#endif
