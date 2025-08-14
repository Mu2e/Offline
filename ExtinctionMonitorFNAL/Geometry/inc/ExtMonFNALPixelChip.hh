//  A pixel chip.
//
// Andrei Gaponenko, 2012

#ifndef EXTMONFNALPIXELCHIP_HH
#define EXTMONFNALPIXELCHIP_HH

#include "canvas/Persistency/Common/Wrapper.h"

namespace mu2e {

  namespace ExtMonFNAL { class ExtMon; }
  class ExtMonFNALModule;

  class ExtMonFNALPixelChip {
  public:

    unsigned int nColumns() const;
    double xPitch() const;
    double xPitch_Edge() const;
    double xPitch_Mid() const;

    unsigned int nRows() const;
    double yPitch() const;

    unsigned int nPixels() const;

  private:
    // Required by genreflex persistency
    ExtMonFNALPixelChip() {}
    template<class T> friend class art::Wrapper;

    friend class ExtMonFNAL::ExtMon;
    friend class ExtMonFNALModule;
  };

} // namespace mu2e

#endif/*EXTMONFNALPIXELCHIP_HH*/
