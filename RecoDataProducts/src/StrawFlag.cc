//
// Class to describe flag bits used for straw materials
//
// Original author David Brown
//
// Mu2e includes
#include "Offline/RecoDataProducts/inc/StrawFlag.hh"
#include <stdexcept>
#include <iostream>
#include <stdio.h>

namespace mu2e {
  unsigned StrawFlagDetail::_maxTrkId(7);

  std::string const& StrawFlagDetail::typeName() {
    static std::string type("StrawFlag");
    return type;
  }

  std::map<std::string,StrawFlagDetail::mask_type> const& StrawFlagDetail::bitNames() {
    static std::map<std::string,mask_type> bitnames;
    if(bitnames.size()==0){
      bitnames[std::string("Active")]         = bit_to_mask(active); // was this material active in the fit?
      bitnames[std::string("HasHit")]         = bit_to_mask(hashit);  // was it associated with a hit?
      bitnames[std::string("ActiveHit")]      = bit_to_mask(activehit); // was its associated hit active?
      bitnames[std::string("DriftHit")]       = bit_to_mask(drifthit); // did its associated hit use drift information in the fit?
    }
    return bitnames;
  }
}
