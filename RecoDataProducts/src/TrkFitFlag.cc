//
// Class to describe flag bits used for straw hits
//
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
      bitnames[std::string("CircleInit")]           = bit_to_mask(circleInit);
      bitnames[std::string("PhiZInit")]           = bit_to_mask(phizInit);
      bitnames[std::string("HitsOK")]             = bit_to_mask(hitsOK);
      bitnames[std::string("CircleOK")]           = bit_to_mask(circleOK);
      bitnames[std::string("PhiZOK")]           = bit_to_mask(phizOK);
      bitnames[std::string("HelixOK")]              = bit_to_mask(helixOK);
      bitnames[std::string("SeedOK")]              = bit_to_mask(seedOK);
      bitnames[std::string("KalmanOK")]              = bit_to_mask(kalmanOK);
      bitnames[std::string("CircleConverged")]           = bit_to_mask(circleConverged);
      bitnames[std::string("PhiZConverged")]           = bit_to_mask(phizConverged);
      bitnames[std::string("HelixConverged")]           = bit_to_mask(helixConverged);
      bitnames[std::string("SeedConverged")]              = bit_to_mask(seedConverged);
      bitnames[std::string("KalmanConverged")]              = bit_to_mask(kalmanConverged);
      bitnames[std::string("KalSeedFit")]              = bit_to_mask(KSF);
      bitnames[std::string("KalFinalFit")]              = bit_to_mask(KFF);
      bitnames[std::string("TrkPatRecHelix")]              = bit_to_mask(TPRHelix);
      bitnames[std::string("CalPatRecHelix")]              = bit_to_mask(CPRHelix);
      bitnames[std::string("Straight")]              = bit_to_mask(Straight);
    }
    return bitnames;
  }

}
