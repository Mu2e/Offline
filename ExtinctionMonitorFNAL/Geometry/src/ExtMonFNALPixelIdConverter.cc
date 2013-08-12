// Andrei Gaponenko, 2012

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPixelIdConverter.hh"

#include <cassert>

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPixelChip.hh"

namespace mu2e {

  //================================================================
  ExtMonFNALPixelIdConverter::ExtMonFNALPixelIdConverter(const mu2e::ExtMonFNAL::ExtMon& extmon)
    : extmon_(&extmon)
    , nPlanes_(extmon.nplanes())
    , module_(extmon.module())
    , chip_(extmon.chip())
    , totalNumberOfPixels_(extmon.nplanes() *
                           extmon.module().nxChips() * extmon.module().nyChips() *
                           extmon.chip().nColumns() * extmon.chip().nRows())
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
    unsigned int globalModuleNumber = (globalChipNumber + 1) / 2;
    ExtMonFNALModuleDenseId did(globalModuleNumber);

    ExtMonFNALModuleIdConverter con(*extmon_);
    ExtMonFNALModuleId mid = con.moduleId(did);

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
