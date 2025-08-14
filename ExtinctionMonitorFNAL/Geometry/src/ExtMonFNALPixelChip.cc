// Andrei Gaponenko, 2012

#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPixelChip.hh"
#include "CLHEP/Units/SystemOfUnits.h"

namespace mu2e {

  // FE-I4 parameters
  unsigned int ExtMonFNALPixelChip::nColumns() const { return 80; }
  double ExtMonFNALPixelChip::xPitch() const { return 250 * CLHEP::micrometer; }
  double ExtMonFNALPixelChip::xPitch_Edge() const { return 500 * CLHEP::micrometer; }
  double ExtMonFNALPixelChip::xPitch_Mid() const  { return 450 * CLHEP::micrometer; }


  unsigned int ExtMonFNALPixelChip::nRows() const { return 336; }
  double ExtMonFNALPixelChip::yPitch() const { return 50 * CLHEP::micrometer; }

  unsigned int ExtMonFNALPixelChip::nPixels() const { return nRows() * nColumns(); }

} // namespace mu2e
