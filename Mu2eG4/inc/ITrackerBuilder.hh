/*
 * ITrackerBuilder.h
 *
 *  Created on: Feb 11, 2010
 *      Author: tassiell
 */

#ifndef Mu2eG4_ITrackerBuilder_hh
#define Mu2eG4_ITrackerBuilder_hh

#include "ITrackerGeom/inc/ITracker.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

// G4 includes
#include "G4LogicalVolume.hh"

namespace mu2e {

  class ITrackerBuilder {
  public:
    static VolumeInfo constructTracker( G4LogicalVolume* mother, double zOff );
  private:
    static VolumeInfo buildWire(float radius, float length, char *shapeName, char *volName, const std::vector<std::string> &materialName, const std::vector<double> &thicknesses, bool activeWireSD=false, bool isSense=false);
    static VolumeInfo buildWall(Wall *wall, ITracker::EnCapType endcapType);
    static double constructSpiderWeb(G4LogicalVolume* localMother, SimpleConfig const& config);
    static void   constructWireAnchoring(G4LogicalVolume* localMother, SimpleConfig const& config, double spdWebBaseExcess=0.0);
    static void   constructSignalCables(G4LogicalVolume* localMother, SimpleConfig const& config);
    static void   constructHvCables(G4LogicalVolume* localMother, SimpleConfig const& config);

  };

} //namespace mu2e

#endif /* Mu2eG4_ITrackerBuilder_hh */
