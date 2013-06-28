// Andrei Gaponenko, 2012

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPixelIdConverter.hh"

#include <cassert>

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPixelChip.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALModule.hh"

namespace mu2e {

  //================================================================
  ExtMonFNALPixelIdConverter::ExtMonFNALPixelIdConverter(unsigned nPlanes,
                                                         const ExtMonFNALModule& module,
                                                         const ExtMonFNALPixelChip& chip)
    : nPlanes_(nPlanes)
    , module_(module)
    , chip_(chip)
    , totalNumberOfPixels_(nPlanes *
                           module.nxChips() * module.nyChips() *
                           chip.nColumns() * chip.nRows())
  {}

  //================================================================
  ExtMonFNALPixelDenseId ExtMonFNALPixelIdConverter::densePixelNumber(const ExtMonFNALPixelId& id) const {

    unsigned res = (
                    id.chip().module().plane() * module_.nxChips() * module_.nyChips()
                    + id.chip().chipRow() * module_.nxChips() + id.chip().chipCol()
                    ) * chip_.nPixels()
      +
      (
       id.row() * chip_.nColumns() + id.col()
       )
      ;

    assert(res < totalNumberOfPixels_);

    return ExtMonFNALPixelDenseId(res);
  }

  //================================================================
  ExtMonFNALPixelId ExtMonFNALPixelIdConverter::pixelId(ExtMonFNALPixelDenseId pix) const {

    assert(pix.number() < totalNumberOfPixels_);

    unsigned int globalChipNumber = pix.number() / chip_.nPixels();

    // WARNING: this code will break if any geometry changes are made.  It assumes 4 planes in the upstream and downstream stack.  
    // All upstream planes have 2 modules and downstream planes have 3 modules.
   
 unsigned int globalModuleNumber = (globalChipNumber + 1) / 2;
    unsigned int globalPlaneNumber = globalModuleNumber <= 8 ? ( (globalModuleNumber + 1) / 2 ) : ( (globalModuleNumber / 3) + 2 );
    
    mu2e::ExtMonFNALModuleId::Rotation rot = globalPlaneNumber <= 4 ? 
      mu2e::ExtMonFNALModuleId::Rotation::UP : ( globalModuleNumber % 3 == 2 ? 
                                                 mu2e::ExtMonFNALModuleId::Rotation::DOWN : mu2e::ExtMonFNALModuleId::Rotation::UP );
   
    mu2e::ExtMonFNALModuleId::Side side = globalPlaneNumber <= 4 ? 
      ( globalModuleNumber % 2 == 1 ? 
        mu2e::ExtMonFNALModuleId::Side::FRONT : mu2e::ExtMonFNALModuleId::Side::BACK ) : ( globalModuleNumber % 3 == 1 ? 
                                                                                           mu2e::ExtMonFNALModuleId::Side::BACK : mu2e::ExtMonFNALModuleId::Side::FRONT );
    
    ExtMonFNALModuleId mid(globalPlaneNumber, side, rot, globalModuleNumber);

    unsigned int chipInModuleNumber = globalChipNumber % (module_.nxChips() * module_.nyChips());
    unsigned int chipRow = chipInModuleNumber / module_.nxChips();
    unsigned int chipCol = chipInModuleNumber % module_.nxChips();

    ExtMonFNALChipId cid(mid, chipCol, chipRow);

    unsigned int pixInChip = pix.number() % chip_.nPixels();
    unsigned int pixRow = pixInChip / chip_.nColumns();
    unsigned int pixCol = pixInChip % chip_.nColumns();

    return ExtMonFNALPixelId(cid, pixCol, pixRow);
  }

} // namespace mu2e
