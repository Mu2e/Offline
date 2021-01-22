//
// Class to describe flag bits used for defining straw (or panel or plane) status
// 
// Original author David Brown (7/2020)
//
// Mu2e includes
#include "DataProducts/inc/StrawStatus.hh"
#include <stdexcept>
#include <iostream>
#include <stdio.h>

namespace mu2e {

  std::string const& StrawStatusDetail::typeName() {
    static std::string type("StrawStatus");
    return type;
  }

  std::map<std::string,StrawStatusDetail::mask_type> const& StrawStatusDetail::bitNames() {
    static std::map<std::string,mask_type> bitnames;
    if(bitnames.size()==0){
      bitnames[std::string("Absent")]         = bit_to_mask(absent);
      bitnames[std::string("NoWire")]         = bit_to_mask(nowire);
      bitnames[std::string("NoHV")]           = bit_to_mask(noHV);
      bitnames[std::string("NoLV")]           = bit_to_mask(noLV);
      bitnames[std::string("NoGas")]          = bit_to_mask(nogas);
      bitnames[std::string("LowGasGain")]     = bit_to_mask(lowgasgain);
      bitnames[std::string("NoPreamp")]      = bit_to_mask(noPreamp);
      bitnames[std::string("NoADC")]	      = bit_to_mask(noADC);
      bitnames[std::string("NoTDC")]          = bit_to_mask(noTDC);
      bitnames[std::string("Sparking")]       = bit_to_mask(sparking);
      bitnames[std::string("Noise")]          = bit_to_mask(noise);
      bitnames[std::string("Pickup")]         = bit_to_mask(pickup);
      bitnames[std::string("Suppress")]       = bit_to_mask(suppress);
    }
    return bitnames;
  }
}
