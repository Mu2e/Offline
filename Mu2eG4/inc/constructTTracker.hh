#ifndef Mu2eG4_constructTTracker_hh
#define Mu2eG4_constructTTracker_hh
//
// Free functions to construct various versions of the TTracker.
//
// $Id: constructTTracker.hh,v 1.7 2013/01/07 04:03:14 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/01/07 04:03:14 $
//
// Original author Rob Kutschke
//

#include "G4Helper/inc/VolumeInfo.hh"

class G4LogicalVolume;

namespace mu2e{

  class SimpleConfig;

  VolumeInfo constructTTrackerv3( VolumeInfo const& mother,
                                  double zOff,
                                  SimpleConfig const& config );

  VolumeInfo constructTTrackerv3Detailed( VolumeInfo const& mother,
                                          double zOff,
                                          SimpleConfig const& config );

}  // end namespace mu2e

#endif /* Mu2eG4_constructTTracker_hh */
