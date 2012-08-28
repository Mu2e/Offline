// Andrei Gaponenko, 2012

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALSensor.hh"

#include <iostream>
#include <cmath>

#include "CLHEP/Units/SystemOfUnits.h"

namespace mu2e {

  // Code in this file (only!) assumes these nubers are (2,2)
  int ExtMonFNALSensor::nxChips() const { return 2; }
  int ExtMonFNALSensor::nyChips() const { return 2; }

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
    const double chipx0 = (icx - 1)*chip_.nColumns()*chip_.xPitch();
    const double chipy0 = (icy - 1)*chip_.nRows()*chip_.yPitch();

    // Pixel column and row numbers are 1-based for FE-I4
    const int ix = 1 + std::floor((xSensor - chipx0)/chip_.xPitch());
    const int iy = 1 + std::floor((ySensor - chipy0)/chip_.yPitch());

    ExtMonFNALPixelId res =
      (0 < ix)&&(ix <= chip_.nColumns())&&(0 < iy)&&(iy <= chip_.nRows()) ?
      ExtMonFNALPixelId(cid, ix, iy) :
      ExtMonFNALPixelId();

    return res;
  }

  //================================================================
  CLHEP::Hep2Vector ExtMonFNALSensor::sensorCoordinates(const ExtMonFNALPixelId& id) const
  {
    // Assume no gaps between chips
    double chipx0 = (id.chip().chipCol() - 1)*chip_.nColumns()*chip_.xPitch();
    double chipy0 = (id.chip().chipRow() - 1)*chip_.nRows()*chip_.yPitch();

    // Subtract 0.5 to get pixel center in the 1-based numbering convention
    const double xSensor = chipx0 + chip_.xPitch()*(id.col() - 0.5);
    const double ySensor = chipy0 + chip_.yPitch()*(id.row() - 0.5);

    return CLHEP::Hep2Vector(xSensor, ySensor);
  }

  //================================================================

} // namespace mu2e
