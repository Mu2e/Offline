// A class to find pixel neighbors
//
// Andrei Gaponenko, 2012

#ifndef ExtinctionMonitorFNAL_Reconstruction_PixelNeighbors_hh
#define ExtinctionMonitorFNAL_Reconstruction_PixelNeighbors_hh

#include <vector>
#include "DataProducts/inc/ExtMonFNALPixelId.hh"

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPixelChip.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALModule.hh"

namespace mu2e {

  class PixelNeighbors {
  public:
    typedef std::vector<ExtMonFNALPixelId>  Collection;

    PixelNeighbors(const ExtMonFNALModule& module, const ExtMonFNALPixelChip& chip);

    Collection neighbors(const ExtMonFNALPixelId& id) const;

  private:
    ExtMonFNALModule module_;
    ExtMonFNALPixelChip chip_;
  };

} // namespace mu2e

#endif/*ExtinctionMonitorFNAL_Reconstruction_PixelNeighbors_hh*/
