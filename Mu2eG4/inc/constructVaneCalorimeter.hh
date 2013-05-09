#ifndef Mu2eG4_constructVaneCalorimeter_hh
#define Mu2eG4_constructVaneCalorimeter_hh
//
// Free function to create the calorimeter.
//
// $Id: constructVaneCalorimeter.hh,v 1.3 2013/05/09 23:14:14 echenard Exp $
// $Author: echenard $
// $Date: 2013/05/09 23:14:14 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) Arguments are:
//    1 - pointer to the mother logical volume.
//    2 - geometry file

// Mu2e includes.
#include "G4Helper/inc/VolumeInfo.hh"

namespace mu2e {

  class SimpleConfig;

  VolumeInfo constructVaneCalorimeter( VolumeInfo const&  mother,SimpleConfig const& config );

}

#endif /* Mu2eG4_constructVaneCalorimeter_hh */
