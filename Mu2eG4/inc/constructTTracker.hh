#ifndef constructTTracker_HH
#define constructTTracker_HH
//
// Free functions to construct various versions of the TTracker.
//
// $Id: constructTTracker.hh,v 1.2 2010/07/19 22:38:44 genser Exp $
// $Author: genser $
// $Date: 2010/07/19 22:38:44 $
//
// Original author Rob Kutschke
//

#include "Mu2eG4/inc/VolumeInfo.hh"

class G4LogicalVolume;

namespace mu2e{

  class SimpleConfig;

  VolumeInfo constructTTrackerv1( G4LogicalVolume* mother,
                                  double zOff,
                                  SimpleConfig const& config );

  VolumeInfo constructTTrackerv2( G4LogicalVolume* mother,
                                  double zOff,
                                  SimpleConfig const& config );


}  // end namespace mu2e

#endif
