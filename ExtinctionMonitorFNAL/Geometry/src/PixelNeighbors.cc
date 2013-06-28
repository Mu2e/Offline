// Andrei Gaponenko, 2012

#include "ExtinctionMonitorFNAL/Geometry/inc/PixelNeighbors.hh"

namespace mu2e {

  //================================================================
  PixelNeighbors::PixelNeighbors(const ExtMonFNALModule& module,
                                 const ExtMonFNALPixelChip& chip)
    : module_(module)
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
    else if(0 < id.chip().chipCol()) {
      ExtMonFNALChipId chip(id.chip().module(), id.chip().chipCol() - 1, id.chip().chipRow());
      res.push_back(ExtMonFNALPixelId(chip, chip_.nColumns()-1, id.row()));
    }

    if(id.col() + 1 < chip_.nColumns()) {
      res.push_back(ExtMonFNALPixelId(id.chip(), id.col()+1, id.row()));
    }
    else if(id.chip().chipCol() + 1 < module_.nxChips()) {
      ExtMonFNALChipId chip(id.chip().module(), id.chip().chipCol() + 1, id.chip().chipRow());
      res.push_back(ExtMonFNALPixelId(chip, 0, id.row()));
    }

    if(0 < id.row()) {
      res.push_back(ExtMonFNALPixelId(id.chip(), id.col(), id.row()-1));
    }
    else if(0 < id.chip().chipRow()) {
      ExtMonFNALChipId chip(id.chip().module(), id.chip().chipCol(), id.chip().chipRow() - 1);
      res.push_back(ExtMonFNALPixelId(chip, id.col(), chip_.nRows()-1));
    }

    if(id.row() + 1 < chip_.nRows()) {
      res.push_back(ExtMonFNALPixelId(id.chip(), id.col(), id.row()+1));
    }
    else if(id.chip().chipRow() + 1 < module_.nyChips()) {
      ExtMonFNALChipId chip(id.chip().module(), id.chip().chipCol(), id.chip().chipRow() + 1);
      res.push_back(ExtMonFNALPixelId(chip, id.col(), 0));
    }

    return res;
  }

  //================================================================

} // namespace mu2e
