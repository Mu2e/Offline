// Andrei Gaponenko, 2012

#include "Mu2eBuildingGeom/inc/BuildingBasicsMaker.hh"

#include <iostream>

#include "cetlib/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "ConfigTools/inc/SimpleConfig.hh"

#include "Mu2eBuildingGeom/inc/BuildingBasics.hh"

namespace mu2e {

  std::unique_ptr<BuildingBasics> BuildingBasicsMaker::make(const SimpleConfig& c)
  {
    std::unique_ptr<BuildingBasics> b(new BuildingBasics());

    b->detectorHallFloorTopY_ = c.getDouble("yOfFloorSurface.below.mu2eOrigin");
    b->detectorHallInsideFullHeight_ = c.getDouble("hall.insideFullHeight");
    b->detectorHallCeilingThickness_ = c.getDouble("hall.ceilingThick");
    b->detectorHallInnerTSCeilingThickness_ = c.getDouble("hall.innerTSCeilingThick");
    b->detectorHallFloorThickness_ = c.getDouble("hall.floorThick");
    b->detectorHallFloorTopDepthBelowGrade_ = c.getDouble("hall.floorTopDepthBelowGrade");

    return b;
  }
}
