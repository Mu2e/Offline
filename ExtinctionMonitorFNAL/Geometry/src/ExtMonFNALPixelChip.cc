// Andrei Gaponenko, 2012

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPixelChip.hh"
#include "CLHEP/Units/SystemOfUnits.h"

namespace mu2e {

  // FE-I4 parameters
  int ExtMonFNALPixelChip::nColumns() const { return 80; }
  double ExtMonFNALPixelChip::xPitch() const { return 250 * CLHEP::micrometer; }

  int ExtMonFNALPixelChip::nRows() const { return 336; }
  double ExtMonFNALPixelChip::yPitch() const { return 50 * CLHEP::micrometer; }

} // namespace mu2e
