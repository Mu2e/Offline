#ifndef Mu2eG4_constructTTracker_hh
#define Mu2eG4_constructTTracker_hh
//
// Free functions to construct various versions of the TTracker.
//
// $Id: constructTTracker.hh,v 1.8 2013/12/20 20:08:21 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/12/20 20:08:21 $
//
// Original author Rob Kutschke
//


class G4LogicalVolume;

namespace mu2e{

  class VOlumeInfo;
  class SimpleConfig;

  VolumeInfo constructTTrackerv3( VolumeInfo const& parent,
                                  SimpleConfig const& config );

  VolumeInfo constructTTrackerv3Detailed( VolumeInfo const& parent,
                                          SimpleConfig const& config );

  VolumeInfo constructTTrackerv5( VolumeInfo const& parent,
                                 SimpleConfig const& config );

}  // end namespace mu2e

#endif /* Mu2eG4_constructTTracker_hh */
