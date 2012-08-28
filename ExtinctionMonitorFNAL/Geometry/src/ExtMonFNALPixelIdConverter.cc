// Andrei Gaponenko, 2012

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPixelIdConverter.hh"

#include <cassert>

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPixelChip.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALSensor.hh"

namespace mu2e {

  //================================================================
  ExtMonFNALPixelIdConverter::ExtMonFNALPixelIdConverter(unsigned nSensorPlanes,
                                                         const ExtMonFNALSensor& sensor,
                                                         const ExtMonFNALPixelChip& chip)
    : nSensorPlanes_(nSensorPlanes)
    , sensor_(sensor)
    , chip_(chip)
    , totalNumberOfPixels_(nSensorPlanes *
                           sensor.nxChips() * sensor.nyChips() *
                           chip.nColumns() * chip.nRows())
  {}

  //================================================================
  ExtMonFNALPixelDenseId ExtMonFNALPixelIdConverter::densePixelNumber(const ExtMonFNALPixelId& id) const {

    unsigned res = (
                    id.chip().sensor().plane() * sensor_.nxChips() * sensor_.nyChips()
                    + id.chip().chipRow() * sensor_.nxChips() + id.chip().chipCol()
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

    ExtMonFNALSensorId sid(globalChipNumber/(sensor_.nxChips() * sensor_.nyChips()));

    unsigned int chipInSensorNumber = globalChipNumber % (sensor_.nxChips() * sensor_.nyChips());
    unsigned int chipRow = chipInSensorNumber / sensor_.nxChips();
    unsigned int chipCol = chipInSensorNumber % sensor_.nxChips();

    ExtMonFNALChipId cid(sid, chipCol, chipRow);

    unsigned int pixInChip = pix.number() % chip_.nPixels();
    unsigned int pixRow = pixInChip / chip_.nColumns();
    unsigned int pixCol = pixInChip % chip_.nColumns();

    return ExtMonFNALPixelId(cid, pixCol, pixRow);
  }

} // namespace mu2e
