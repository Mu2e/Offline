// Evan Schiewe, 2013

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPlane.hh"

#include <iostream>
#include <cmath>

#include "CLHEP/Units/SystemOfUnits.h"

namespace mu2e {
  
  CLHEP::Hep3Vector ExtMonFNALPlane::planeCoordinates(const ExtMonFNALPixelId& id) const
  {

    CLHEP::Hep2Vector pixelCordinatesInModule = this->module().moduleCoordinates(id);
 
   // checking for module rotation; DOWN is 180 degrees rotated from UP, so coordinates are flipped
    const int rotz = this->module_rotation()[(id.chip().module().number())] == 0 ? 1 : -1;

    // check for module rotation about y
    const int roty = this->module_zoffset()[(id.chip().module().number())] > 0 ? 1 : -1;

    const double xModule = this->module_xoffset()[(id.chip().module().number())] + pixelCordinatesInModule.x()*rotz*roty;
    const double yModule = this->module_yoffset()[(id.chip().module().number())] + pixelCordinatesInModule.y()*rotz;
    // zModule in the geometry file is +1 or -1, representing position on the front or back of the plane.  It is not
    // the actual offset, unlike xModule and yModule. Here we convert it to actual coordinates.
    const double zModule = this->module_zoffset()[(id.chip().module().number())] 
	* (this->module().sensorHalfSize()[2] + 2*this->module().chipHalfSize()[2] + this->halfSize()[2]);

    return CLHEP::Hep3Vector(xModule, yModule, zModule);
  }

 } // namespace mu2e
 
