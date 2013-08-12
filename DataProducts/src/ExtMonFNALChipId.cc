#include "DataProducts/inc/ExtMonFNALChipId.hh"

namespace mu2e {
  ExtMonFNALChipId::ExtMonFNALChipId(const ExtMonFNALModuleId& module,
                                     unsigned int chipCol,
                                     unsigned int chipRow)
    : module_(module)
    , chipCol_(chipCol)
    , chipRow_(chipRow)
  {}


  std::ostream& operator<<( std::ostream& os, const ExtMonFNALChipId& id) {
    return os<<"ExtMonFNALChipId("<<id.module().plane()
             <<","<<id.module().number()
             <<","<<id.chipCol()
             <<","<<id.chipRow()
             <<" )";
  }
}
