#include "Offline/Mu2eKinKal/inc/WHSMask.hh"
#include <stdexcept>
#include <iostream>
#include <stdio.h>

namespace mu2e {

  std::string const& WHSMaskDetail::typeName() {
    static std::string type("WHSMask");
    return type;
  }

  std::map<std::string,WHSMaskDetail::mask_type> const& WHSMaskDetail::bitNames() {
    static std::map<std::string,mask_type> bitnames;
    if(bitnames.size()==0){
      bitnames[std::string("Inactive")]           = bit_to_mask(inactive);
      bitnames[std::string("Null")]           = bit_to_mask(null);
      bitnames[std::string("Drift")]             = bit_to_mask(drift);
      }
    return bitnames;
  }

}
