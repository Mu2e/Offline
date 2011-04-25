#ifndef constructCalorimeter_HH
#define constructCalorimeter_HH
//
// Free function to create the calorimeter.
//
// $Id: constructCalorimeter.hh,v 1.4 2011/04/25 19:16:49 genser Exp $
// $Author: genser $
// $Date: 2011/04/25 19:16:49 $
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

#endif
