// Evan Schiewe, 2013

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALModule.hh"

#include <iostream>
#include <cmath>

#include "CLHEP/Units/SystemOfUnits.h"

namespace mu2e {

  // Code in this file (only!) assumes these numbers are (2,1) 
  unsigned int ExtMonFNALModule::nxChips() const { return 2; }
  unsigned int ExtMonFNALModule::nyChips() const { return 1; }

  //================================================================
  ExtMonFNALPixelId ExtMonFNALModule::findPixel(ExtMonFNALModuleId mid,
                                                double xSensor,
                                                double ySensor) const
  {
    // We assume there are 2x1 chips per module
    const unsigned icx = xSensor < 0 ? 0 : 1;
    const unsigned icy = 0;
    ExtMonFNALChipId cid(mid, icx, icy);

    // Assume no gaps between chips
    const double chipx0 = (icx - 1.)*chip_.nColumns()*chip_.xPitch();
    const double chipy0 = (icy - 1.)*chip_.nRows()*chip_.yPitch();    

    
    // enum Rotation {UP, DOWN};
    // const int rot = mid.rot() == UP ? 1 : -1;

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
  CLHEP::Hep2Vector ExtMonFNALModule::moduleCoordinates(const ExtMonFNALPixelId& id) const
  {
    // Assume no gaps between chips
    double chipx0 = (id.chip().chipCol() - 1.)*chip_.nColumns()*chip_.xPitch();
    double chipy0 = (id.chip().chipRow() - 1.)*chip_.nRows()*chip_.yPitch();

    // Add 0.5 to get pixel center in the 0-based numbering convention
    const double xModule = chipx0 + chip_.xPitch()*(id.col() + 0.5);
    const double yModule = chipy0 + chip_.yPitch()*(id.row() + 0.5);

    return CLHEP::Hep2Vector(xModule, yModule);
  }

  //================================================================

} // namespace mu2e
