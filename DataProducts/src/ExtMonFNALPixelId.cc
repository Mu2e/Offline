#include "DataProducts/inc/ExtMonFNALPixelId.hh"

namespace mu2e {
  ExtMonFNALPixelId::ExtMonFNALPixelId(const ExtMonFNALChipId& chip,
                                       int col,
                                       int row)
    : chip_(chip)
    , col_(col)
    , row_(row)
  {}

  std::ostream& operator<<( std::ostream& os, const ExtMonFNALPixelId& id) {
    return os<<"ExtMonFNALPixelId("<<id.chip().sensor().plane()
             <<", "<<id.chip().chipCol()
             <<", "<<id.chip().chipRow()
             <<", "<<id.col()
             <<", "<<id.row()
             <<" )";
  }
}
