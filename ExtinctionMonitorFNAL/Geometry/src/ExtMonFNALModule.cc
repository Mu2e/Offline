// Evan Schiewe, 2013

#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALModule.hh"

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
    ExtMonFNALPixelId res; // default constructed - invalid hit

    const double chipXPitch = (chip_.nColumns()-2)*chip_.xPitch() + chip_.xPitch_Edge() + chip_.xPitch_Mid();
    const double chipYPitch = chip_.nRows()*chip_.yPitch();

    const int icx = std::floor(xSensor/chipXPitch + nxChips()/2.);
    const int icy = std::floor(ySensor/chipYPitch + nyChips()/2.);

    if((0 <= icx)&&(unsigned(icx) < nxChips())&&(0 <= icy)&&(unsigned(icy) < nyChips())) {

        ExtMonFNALChipId cid(mid, icx, icy);

        // Assume no gaps between chips
        // x0 and y0 are the coordinates of the bottom left chip corner in the module frame
        const double chipx0 = (icx - nxChips()/2.)*chip_.nColumns()*chip_.xPitch();
        const double chipy0 = (icy - nyChips()/2.)*chip_.nRows()*chip_.yPitch();

        // Zero based pixel column and row numbers for the offline identifier
        //const int ix = std::floor((xSensor - chipx0)/chip_.xPitch());
        const int iy = std::floor((ySensor - chipy0)/chip_.yPitch());
        int ix_=-1;

        if(xSensor > (chip_.nColumns()-2)*chip_.xPitch()+chip_.xPitch_Mid())
       {
          ix_ = chip_.nColumns() - 1;
       }

        else if(xSensor>chip_.xPitch_Mid() && xSensor<= (chip_.nColumns()-2)*chip_.xPitch()+chip_.xPitch_Mid())
       {
          ix_ = floor((xSensor - chipx0 - chip_.xPitch_Mid() + chip_.xPitch())/chip_.xPitch());
       }

        else if(xSensor>=0 && xSensor<=chip_.xPitch_Mid())
       {
          ix_ = 0;
       }

        else if(xSensor<0 && xSensor>=-chip_.xPitch_Mid())
       {
          ix_ = chip_.nColumns() - 1;
       }

        else if(xSensor<-chip_.xPitch_Mid() && xSensor>=  -(chip_.nColumns()-2)*chip_.xPitch()-chip_.xPitch_Mid() )
       {
          ix_ = floor((xSensor - chipx0)/chip_.xPitch_Mid()) - 1;
       }

        else if (xSensor < -(chip_.nColumns()-2)*chip_.xPitch()-chip_.xPitch_Mid())
       {
          ix_ = 0;
       }

        const int ix = ix_;

        res = (0 <= ix)&&(unsigned(ix) < chip_.nColumns())&&(0 <= iy)&&(unsigned(iy) < chip_.nRows()) ?
           ExtMonFNALPixelId(cid, ix, iy) :
           ExtMonFNALPixelId();
    }

    return res;
  }

  //================================================================
  CLHEP::Hep2Vector ExtMonFNALModule::moduleCoordinates(const ExtMonFNALPixelId& id) const
  {
    // Assume no gaps between chips
    double chipx0 = (id.chip().chipCol() - nxChips()/2.)*chip_.nColumns()*chip_.xPitch();
    double chipy0 = (id.chip().chipRow() - nyChips()/2.)*chip_.nRows()*chip_.yPitch();

    // Add 0.5 to get pixel center in the 0-based numbering convention
    const double xModule = chipx0 + chip_.xPitch()*(id.col() + 0.5);
    const double yModule = chipy0 + chip_.yPitch()*(id.row() + 0.5);

    return CLHEP::Hep2Vector(xModule, yModule);
  }

  //================================================================

} // namespace mu2e
