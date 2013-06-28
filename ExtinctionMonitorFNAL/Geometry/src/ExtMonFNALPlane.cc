// Evan Schiewe, 2013

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPlane.hh"

#include <iostream>
#include <cmath>

#include "CLHEP/Units/SystemOfUnits.h"

namespace mu2e {
  
  CLHEP::Hep3Vector ExtMonFNALPlane::planeCoordinates(const ExtMonFNALPixelId& id) const
  {
    // Assume no gaps between chips
    double chipx0 = (id.chip().chipCol() - 1.)*80*250*CLHEP::micrometer;
    double chipy0 = (id.chip().chipRow() - 1.)*336*50*CLHEP::micrometer;

    // checking for module rotation; DOWN is 180 degrees rotated from UP, so coordinates are flipped
    enum Rotation {UP, DOWN};
    const int rot = id.chip().module().rot() == UP ? 1 : -1;

    // Add 0.5 to get pixel center in the 0-based numbering convention
    const double xModule = (chipx0 + 250*CLHEP::micrometer*(id.col() + 0.5)) * rot;
    const double yModule = (chipy0 + 50*CLHEP::micrometer*(id.row() + 0.5)) * rot;

    enum Side {FRONT, BACK};
    const double zModule = (id.chip().module().side() == FRONT ? .3 : -.3);

    return CLHEP::Hep3Vector(xModule, yModule, zModule);
  }

 } // namespace mu2e
 
