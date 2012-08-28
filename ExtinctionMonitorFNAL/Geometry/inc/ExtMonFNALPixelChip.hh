//  A pixel chip.
//
// Andrei Gaponenko, 2012

#ifndef EXTMONFNALPIXELCHIP_HH
#define EXTMONFNALPIXELCHIP_HH

#include "art/Persistency/Common/Wrapper.h"

namespace mu2e {

  namespace ExtMonFNAL { class ExtMon; }
  class ExtMonFNALSensor;

  class ExtMonFNALPixelChip {
  public:

    int    nColumns() const;
    double xPitch() const;

    int    nRows() const;
    double yPitch() const;

  private:
    // Required by genreflex persistency
    ExtMonFNALPixelChip() {}
    template<class T> friend class art::Wrapper;

    friend class ExtMonFNAL::ExtMon;
    friend class ExtMonFNALSensor;
  };

} // namespace mu2e

#endif/*EXTMONFNALPIXELCHIP_HH*/
