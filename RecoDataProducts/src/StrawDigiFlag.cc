//
// Class to describe flag bits for StrawDigis
//
//
// Original author David Brown
//
// Mu2e includes
#include "Offline/RecoDataProducts/inc/StrawDigiFlag.hh"
#include <stdexcept>
#include <iostream>
#include <stdio.h>

namespace mu2e {
  std::string const& StrawDigiFlagDetail::typeName() {
    static std::string type("StrawDigiFlag");
    return type;
  }

  std::map<std::string,StrawDigiFlagDetail::mask_type> const& StrawDigiFlagDetail::bitNames() {
    static std::map<std::string,mask_type> bitnames;
    if(bitnames.size()==0){
      bitnames[std::string("EnergySelection")]      = bit_to_mask(energysel);
      bitnames[std::string("dE/dxSelection")]      = bit_to_mask(dedxsel);
      bitnames[std::string("processed")]      = bit_to_mask(dedxsel);
    }
    return bitnames;
  }
}
