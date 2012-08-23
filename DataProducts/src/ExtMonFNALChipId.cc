#include "DataProducts/inc/ExtMonFNALChipId.hh"

namespace mu2e {
  ExtMonFNALChipId::ExtMonFNALChipId(const ExtMonFNALSensorId& sensor,
                                     int chipCol,
                                     int chipRow)
    : sensor_(sensor)
    , chipCol_(chipCol)
    , chipRow_(chipRow)
  {}


  std::ostream& operator<<( std::ostream& os, const ExtMonFNALChipId& id) {
    return os<<"ExtMonFNALChipId("<<id.sensor().plane()
             <<","<<id.chipCol()
             <<","<<id.chipRow()
             <<" )";
  }
}
