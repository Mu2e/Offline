#ifndef constructLTracker_HH
#define constructLTracker_HH
//
// Free functions to construct various versions of the LTracker.
//
// $Id: constructLTracker.hh,v 1.2 2010/12/22 17:37:57 genser Exp $
// $Author: genser $
// $Date: 2010/12/22 17:37:57 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) The versions are:
//    0 - just a mother volume; no straws (for graphics only).
//    1 - straws with no substructure
//    2 - straws in assembly volumes per sector
//    3 - separate physical mother volumes per sector.

#include "G4Helper/inc/VolumeInfo.hh"

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
