#ifndef Mu2eG4_constructVaneCalorimeter_hh
#define Mu2eG4_constructVaneCalorimeter_hh
//
// Free function to create the calorimeter.
//
// $Id: constructVaneCalorimeter.hh,v 1.2 2013/03/05 20:33:25 aluca Exp $
// $Author: aluca $
// $Date: 2013/03/05 20:33:25 $
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

  VolumeInfo constructVaneCalorimeter( VolumeInfo const&   mother,
                                 double              zOffset,
                                 SimpleConfig const& config );

}

#endif /* Mu2eG4_constructVaneCalorimeter_hh */
