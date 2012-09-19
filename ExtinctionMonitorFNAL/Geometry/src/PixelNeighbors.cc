// Andrei Gaponenko, 2012

#include "ExtinctionMonitorFNAL/Geometry/inc/PixelNeighbors.hh"

namespace mu2e {

  //================================================================
  PixelNeighbors::PixelNeighbors(const ExtMonFNALSensor& sensor,
                                 const ExtMonFNALPixelChip& chip)
    : sensor_(sensor)
    , chip_(chip)
  {}

  //================================================================
  PixelNeighbors::Collection
  PixelNeighbors::neighbors(const ExtMonFNALPixelId& id) const {
    Collection res;

    // This is the definition of a neighbor pixel.  Neighboring pixels
    // have a common side, not just touch in the corner.  For now
    // neighbors should be in the same chip.  We may want to extend
    // the definition to accept hits on adjacent chips.

    if(0 < id.col()) {
      res.push_back(ExtMonFNALPixelId(id.chip(), id.col()-1, id.row()));
    }

    if(id.col() + 1 < chip_.nColumns()) {
      res.push_back(ExtMonFNALPixelId(id.chip(), id.col()+1, id.row()));
    }

    if(0 < id.row()) {
      res.push_back(ExtMonFNALPixelId(id.chip(), id.col(), id.row()-1));
    }

    if(id.row() + 1 < chip_.nRows()) {
      res.push_back(ExtMonFNALPixelId(id.chip(), id.col(), id.row()+1));
    }

    return res;
  }

  //================================================================

} // namespace mu2e
