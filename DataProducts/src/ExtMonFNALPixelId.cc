#include "DataProducts/inc/ExtMonFNALPixelId.hh"

namespace mu2e {
  ExtMonFNALPixelId::ExtMonFNALPixelId(const ExtMonFNALChipId& chip,
                                       unsigned int col,
                                       unsigned int row)
    : chip_(chip)
    , col_(col)
    , row_(row)
  {}

  std::ostream& operator<<( std::ostream& os, const ExtMonFNALPixelId& id) {
    return os<<"ExtMonFNALPixelId("<<id.chip().module()
             <<", "<<id.chip().chipCol()
             <<", "<<id.chip().chipRow()
             <<", "<<id.col()
             <<", "<<id.row()
             <<" )";
  }
}
