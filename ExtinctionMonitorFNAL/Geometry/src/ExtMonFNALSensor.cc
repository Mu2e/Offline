// Andrei Gaponenko, 2012

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALSensor.hh"

#include <iostream>
#include <cmath>

#include "CLHEP/Units/SystemOfUnits.h"

namespace mu2e {
  namespace {
    // FE-I4 parameters
    const double xPitch   = 250 * CLHEP::micrometer;
    const int    nColumns = 80;
    const double chipdx = nColumns * xPitch;

    const double yPitch =   50 * CLHEP::micrometer;
    const int    nRows  =  336;
    const double chipdy = nRows * yPitch;
  }

  //================================================================
  ExtMonFNALPixelId ExtMonFNALSensor::findPixel(ExtMonFNALSensorId sid,
                                                double xSensor,
                                                double ySensor) const
  {
    // We assume there are 2x2 chips per sensor, one per coordinate quadrant.
    const int icx = xSensor < 0 ? 0 : 1;
    const int icy = ySensor < 0 ? 0 : 1;
    ExtMonFNALChipId cid(sid, icx, icy);

    // Assume no gaps between chips
    const double chipx0 = (icx - 1)*chipdx;
    const double chipy0 = (icy - 1)*chipdy;

    // Pixel column and row numbers are 1-based for FE-I4
    const int ix = 1 + std::floor((xSensor - chipx0)/xPitch);
    const int iy = 1 + std::floor((ySensor - chipy0)/yPitch);

    ExtMonFNALPixelId res =
      (0 < ix)&&(ix <= nColumns)&&(0 < iy)&&(iy <= nRows) ?
      ExtMonFNALPixelId(cid, ix, iy) :
      ExtMonFNALPixelId();

    return res;
  }

  //================================================================
  CLHEP::Hep2Vector ExtMonFNALSensor::sensorCoordinates(const ExtMonFNALPixelId& id) const
  {
    // Assume no gaps between chips
    double chipx0 = (id.chip().chipCol() - 1)*chipdx;
    double chipy0 = (id.chip().chipRow() - 1)*chipdy;

    // Subtract 0.5 to get pixel center in the 1-based numbering convention
    const double xSensor = chipx0 + xPitch*(id.col() - 0.5);
    const double ySensor = chipy0 + yPitch*(id.row() - 0.5);

    return CLHEP::Hep2Vector(xSensor, ySensor);
  }

  //================================================================

} // namespace mu2e
