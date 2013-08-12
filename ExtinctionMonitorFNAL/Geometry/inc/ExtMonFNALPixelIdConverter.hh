//  Two-way conversion between pixel Id and sequential ("dense") pixel number.
//
// Andrei Gaponenko, 2012

#ifndef EXTMONFNALPIXELIDCONVERTER_HH
#define EXTMONFNALPIXELIDCONVERTER_HH

#include "DataProducts/inc/ExtMonFNALPixelId.hh"
#include "DataProducts/inc/ExtMonFNALPixelDenseId.hh"
#include "DataProducts/inc/ExtMonFNALModuleDenseId.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPixelChip.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALModule.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALModuleIdConverter.hh"

namespace mu2e {

  class ExtMonFNALPixelIdConverter {
  public:

    // should never be called with invalid id
    ExtMonFNALPixelDenseId densePixelNumber(const ExtMonFNALPixelId& id) const;

    // should never be called with invalid pix
    ExtMonFNALPixelId pixelId(ExtMonFNALPixelDenseId pix) const;

    unsigned int totalNumberOfPixels() const { return totalNumberOfPixels_;  }

    explicit  ExtMonFNALPixelIdConverter(const mu2e::ExtMonFNAL::ExtMon& extmon);

  private:
    const mu2e::ExtMonFNAL::ExtMon *extmon_;
    unsigned nPlanes_;
    ExtMonFNALModule module_;
    ExtMonFNALPixelChip chip_;
    unsigned totalNumberOfPixels_;
  };

} // namespace mu2e

#endif/*EXTMONFNALPIXELIDCONVERTER_HH*/
