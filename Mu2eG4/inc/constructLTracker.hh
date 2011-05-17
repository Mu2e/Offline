#ifndef Mu2eG4_constructLTracker_hh
#define Mu2eG4_constructLTracker_hh
//
// Free functions to construct various versions of the LTracker.
//
// $Id: constructLTracker.hh,v 1.3 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
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

#endif /* Mu2eG4_constructLTracker_hh */
