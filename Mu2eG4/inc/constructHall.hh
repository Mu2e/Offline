#ifndef Mu2eG4_constructHall_hh
#define Mu2eG4_constructHall_hh
//
// Free function to create the hall walls and hall interior inside the earthen overburden.
//
// $Id: constructHall.hh,v 1.4 2012/04/17 19:56:56 gandr Exp $
// $Author: gandr $
// $Date: 2012/04/17 19:56:56 $
//
// Original author: Kyle Knoepfel

// C++ includes
#include <map>
#include <string>

// CLHEP includes
#include "CLHEP/Vector/Rotation.h"

namespace mu2e {

  class ExtrudedSolid;
  class SimpleConfig;
  class VolumeInfo;
  class NotchManager;

  VolumeInfo constructHall(const VolumeInfo& worldInfo, const SimpleConfig& config);

  void constructSolids( const SimpleConfig& config,
                        const VolumeInfo& hallInfo, 
		        const std::map<std::string,ExtrudedSolid>& solidMap,
			const CLHEP::HepRotation& rot,
			const NotchManager& notchMgr);

}

#endif /* Mu2eG4_constructHall_hh */
