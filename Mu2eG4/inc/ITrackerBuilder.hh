/*
 * ITrackerBuilder.h
 *
 *  Created on: Feb 11, 2010
 *      Author: tassiell
 */

#ifndef ITRACKERBUILDER_H_
#define ITRACKERBUILDER_H_

#include "ITrackerGeom/inc/ITracker.hh"
#include "G4Helper/inc/VolumeInfo.hh"

// G4 includes
#include "G4LogicalVolume.hh"

namespace mu2e {

  class ITrackerBuilder {
  public:
    static VolumeInfo constructTracker( G4LogicalVolume* mother, double zOff );
  private:
    static VolumeInfo buildWire(float radius, float length, char *shapeName, char *volName, const std::vector<std::string> &materialName, const std::vector<double> &thicknesses);
    static VolumeInfo buildWall(Wall *wall, ITracker::EnCapType endcapType);

  };

} //namespace mu2e

#endif /* ITRACKERBUILDER_H_ */
