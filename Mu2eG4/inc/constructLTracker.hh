#ifndef constructLTracker_HH
#define constructLTracker_HH
//
// Free functions to construct various versions of the LTracker.
//
// $Id: constructLTracker.hh,v 1.1 2010/04/15 23:01:39 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/15 23:01:39 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) The versions are:
//    0 - just a mother volume; no straws (for graphics only).
//    1 - straws with no substructure
//    2 - straws in assembly volumes per sector
//    3 - separate physical mother volumes per sector.

#include "Mu2eG4/inc/VolumeInfo.hh"

class G4LogicalVolume;

namespace mu2e{

  class SimpleConfig;

  VolumeInfo constructLTrackerv1( G4LogicalVolume* mother,
                                  double zOff,
                                  SimpleConfig const& config );

  VolumeInfo constructLTrackerv2( G4LogicalVolume* mother,
                                  double zOff,
                                  SimpleConfig const& config );

  VolumeInfo constructLTrackerv3( G4LogicalVolume* mother,
                                  double zOff,
                                  SimpleConfig const& config );

}  // end namespace mu2e

#endif
