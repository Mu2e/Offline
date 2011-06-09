#ifndef Mu2eG4_constructTTracker_hh
#define Mu2eG4_constructTTracker_hh
//
// Free functions to construct various versions of the TTracker.
//
// $Id: constructTTracker.hh,v 1.6 2011/06/09 19:58:14 genser Exp $
// $Author: genser $
// $Date: 2011/06/09 19:58:14 $
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

  VolumeInfo constructTTrackerv3( VolumeInfo const& mother,
                                  double zOff,
                                  SimpleConfig const& config );

}  // end namespace mu2e

#endif /* Mu2eG4_constructTTracker_hh */
