#ifndef Mu2eG4_constructHall_hh
#define Mu2eG4_constructHall_hh
//
// Free function to create the hall walls and hall interior inside the earthen overburden.
//
//
// Original author: Kyle Knoepfel

// C++ includes
#include <map>
#include <string>

// CLHEP includes
#include "CLHEP/Vector/Rotation.h"

namespace mu2e {

  class ExtrudedSolid;
  class RotExtrudedSolid;
  class GenericTrap;
  class SimpleConfig;
  class VolumeInfo;
  class NotchHoleManager;

  VolumeInfo constructHall(const VolumeInfo& worldInfo, const SimpleConfig& config);

  void constructSolids( const SimpleConfig& config,
                        const VolumeInfo& hallInfo,
                        const std::map<std::string,ExtrudedSolid>& solidMap,
                        const CLHEP::HepRotation& rot,
                        const NotchHoleManager& notchMgr);

  void constructRotSolids( const SimpleConfig& config,
                        const VolumeInfo& hallInfo,
                        const std::map<std::string,RotExtrudedSolid>& solidMap,
                        const CLHEP::HepRotation& rot,
                        const NotchHoleManager& notchMgr);

  void constructTrapSolids( const SimpleConfig& config,
                            const VolumeInfo& hallInfo,
                            const std::map<std::string,GenericTrap>& solidMap,
                            const CLHEP::HepRotation& rot,
                            const NotchHoleManager& notchMgr);

}

#endif /* Mu2eG4_constructHall_hh */
