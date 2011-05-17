#ifndef Mu2eG4_constructTTracker_hh
#define Mu2eG4_constructTTracker_hh
//
// Free functions to construct various versions of the TTracker.
//
// $Id: constructTTracker.hh,v 1.5 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
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

#endif /* Mu2eG4_constructTTracker_hh */
