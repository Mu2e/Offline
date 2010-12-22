#ifndef constructCalorimeter_HH
#define constructCalorimeter_HH
//
// Free function to create the calorimeter.
//
// $Id: constructCalorimeter.hh,v 1.3 2010/12/22 17:37:57 genser Exp $
// $Author: genser $
// $Date: 2010/12/22 17:37:57 $
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

// Forward references.
class G4LogicalVolume;

namespace mu2e {

  class SimpleConfig;

  void constructCalorimeter( G4LogicalVolume    *mother, 
			     double              zOffset,
			     SimpleConfig const& config );

}

#endif
