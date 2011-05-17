#ifndef Mu2eG4_constructCalorimeter_hh
#define Mu2eG4_constructCalorimeter_hh
//
// Free function to create the calorimeter.
//
// $Id: constructCalorimeter.hh,v 1.5 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
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

  void constructCalorimeter( VolumeInfo const&   mother, 
			     double              zOffset,
			     SimpleConfig const& config );

}

#endif /* Mu2eG4_constructCalorimeter_hh */
