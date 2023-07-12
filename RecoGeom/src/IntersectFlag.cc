#include "Offline/RecoGeom/inc/IntersectFlag.hh"
#include <stdexcept>
#include <iostream>
#include <stdio.h>

namespace mu2e {
  namespace RecoGeom {
    unsigned IntersectFlagDetail::_maxTrkId(7);

    std::string const& IntersectFlagDetail::typeName() {
      static std::string type("IntersectFlag");
      return type;
    }

    std::map<std::string,IntersectFlagDetail::mask_type> const& IntersectFlagDetail::bitNames() {
      static std::map<std::string,mask_type> bitnames;
      if(bitnames.size()==0){
        bitnames[std::string("OnSurface")]               = bit_to_mask(onsurface);
        bitnames[std::string("InBounds")]           = bit_to_mask(inbounds);
      }
      return bitnames;
    }
  }
}
