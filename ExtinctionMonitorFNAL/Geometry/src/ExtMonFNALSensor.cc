// Andrei Gaponenko, 2012

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALSensor.hh"

#include <iostream>
#include <cmath>

#include "CLHEP/Units/SystemOfUnits.h"

namespace mu2e {

  // Code in this file (only!) assumes these nubers are (2,2)
  unsigned int ExtMonFNALSensor::nxChips() const { return 2; }
  unsigned int ExtMonFNALSensor::nyChips() const { return 2; }

  //================================================================
  ExtMonFNALPixelId ExtMonFNALSensor::findPixel(ExtMonFNALSensorId sid,
                                                double xSensor,
                                                double ySensor) const
  {
    // We assume there are 2x2 chips per sensor, one per coordinate quadrant.
    const unsigned icx = xSensor < 0 ? 0 : 1;
    const unsigned icy = ySensor < 0 ? 0 : 1;
    ExtMonFNALChipId cid(sid, icx, icy);

    // Assume no gaps between chips
    const double chipx0 = (icx - 1)*chip_.nColumns()*chip_.xPitch();
    const double chipy0 = (icy - 1)*chip_.nRows()*chip_.yPitch();

    // Zero based pixel column and row numbers for the offline identifier
    const int ix = std::floor((xSensor - chipx0)/chip_.xPitch());
    const int iy = std::floor((ySensor - chipy0)/chip_.yPitch());

    ExtMonFNALPixelId res =
      (0 <= ix)&&(unsigned(ix) < chip_.nColumns())&&(0 <= iy)&&(unsigned(iy) < chip_.nRows()) ?
      ExtMonFNALPixelId(cid, ix, iy) :
      ExtMonFNALPixelId();

    return res;
  }

  //================================================================
  CLHEP::Hep2Vector ExtMonFNALSensor::sensorCoordinates(const ExtMonFNALPixelId& id) const
  {
    // Assume no gaps between chips
    double chipx0 = id.chip().chipCol()*chip_.nColumns()*chip_.xPitch();
    double chipy0 = id.chip().chipRow()*chip_.nRows()*chip_.yPitch();

    // Add 0.5 to get pixel center in the 0-based numbering convention
    const double xSensor = chipx0 + chip_.xPitch()*(id.col() + 0.5);
    const double ySensor = chipy0 + chip_.yPitch()*(id.row() + 0.5);

    return CLHEP::Hep2Vector(xSensor, ySensor);
  }

  //================================================================

} // namespace mu2e
