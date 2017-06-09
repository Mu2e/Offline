//
// Class to describe flag bits used for a background cluster
//
// $Id: BkgClusterFlag.cc,v 1.4 2013/04/04 01:08:20 brownd Exp $
// $Author: brownd $
// $Date: 2013/04/04 01:08:20 $
//
// Original author David Brown
//
// Mu2e includes
#include "RecoDataProducts/inc/BkgClusterFlag.hh"
#include <stdexcept>
#include <iostream>
#include <stdio.h>

namespace mu2e {

  std::string const& BkgClusterFlagDetail::typeName() {
    static std::string type("BkgClusterFlag");
    return type;
  }

  std::map<std::string,BkgClusterFlagDetail::mask_type> const& BkgClusterFlagDetail::bitNames() {
    static std::map<std::string,mask_type> bitnames;
    if(bitnames.size()==0){
      bitnames[std::string("stereo")]           = bit_to_mask(stereo);
      bitnames[std::string("tight")]           = bit_to_mask(tight);
      bitnames[std::string("loose")]           = bit_to_mask(loose);
      bitnames[std::string("refined")]           = bit_to_mask(refined);
      bitnames[std::string("background")]       = bit_to_mask(bkg);
      bitnames[std::string("isolated")]       = bit_to_mask(iso);
   }
    return bitnames;
  }

}
