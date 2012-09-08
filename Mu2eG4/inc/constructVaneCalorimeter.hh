#ifndef Mu2eG4_constructVaneCalorimeter_hh
#define Mu2eG4_constructVaneCalorimeter_hh
//
// Free function to create the calorimeter.
//
// $Id: constructVaneCalorimeter.hh,v 1.1 2012/09/08 02:24:25 echenard Exp $
// $Author: echenard $
// $Date: 2012/09/08 02:24:25 $
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

  void constructVaneCalorimeter( VolumeInfo const&   mother,
                                 double              zOffset,
                                 SimpleConfig const& config );

}

#endif /* Mu2eG4_constructVaneCalorimeter_hh */
