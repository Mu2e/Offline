//
// Class to describe flag bits used for straw hits
//
// $Id: TrkFitFlag.cc,v 1.4 2013/04/04 01:08:20 brownd Exp $
// $Author: brownd $
// $Date: 2013/04/04 01:08:20 $
//
// Original author David Brown
//
// Mu2e includes
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include <stdexcept>
#include <iostream>
#include <stdio.h>

namespace mu2e {

  std::string const& TrkFitFlagDetail::typeName() {
    static std::string type("TrkFitFlag");
    return type;
  }

  std::map<std::string,TrkFitFlagDetail::mask_type> const& TrkFitFlagDetail::bitNames() {
    static std::map<std::string,mask_type> bitnames;
    if(bitnames.size()==0){
      bitnames[std::string("Initialized")]        = bit_to_mask(initialized);
      bitnames[std::string("HitsOK")]             = bit_to_mask(hitsOK);
      bitnames[std::string("CenterOK")]           = bit_to_mask(centerOK);
      bitnames[std::string("RadiusOK")]           = bit_to_mask(radiusOK);
      bitnames[std::string("PhiZOK")]             = bit_to_mask(phizOK);
      bitnames[std::string("FitOK")]              = bit_to_mask(fitOK);
    }
    return bitnames;
  }

}
