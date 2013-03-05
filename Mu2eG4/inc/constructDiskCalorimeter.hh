#ifndef Mu2eG4_constructDiskCalorimeter_hh
#define Mu2eG4_constructDiskCalorimeter_hh
//
// Free function to create the disk calorimeter.
//
//
// Original author Rob Kutschke
//
// Notes:
// 1) Arguments are:
//    1 - pointer to the mother logical volume.
//    2 - z postition of the origin of the Mu2e coordintate system in the
//        frame of the mother.
//    3 - geometry file

// Mu2e includes.
#include "G4Helper/inc/VolumeInfo.hh"

namespace mu2e {

  class SimpleConfig;

  VolumeInfo constructDiskCalorimeter( VolumeInfo const&   mother,
                             double              zOffset,
                             SimpleConfig const& config );

}

#endif /* Mu2eG4_constructDiskCalorimeter_hh */
