//
// Original author David Brown
//
// Mu2e includes
#include "Offline/Mu2eKinKal/inc/KKSHFlag.hh"
#include <stdexcept>
#include <iostream>
#include <stdio.h>

namespace mu2e {
  std::string const& KKSHFlagDetail::typeName() {
    static std::string type("KKSHFlag");
    return type;
  }

  std::map<std::string,KKSHFlagDetail::mask_type> const& KKSHFlagDetail::bitNames() {
    static std::map<std::string,mask_type> bitnames;
    if(bitnames.size()==0){
      bitnames[std::string("TOT")]              = bit_to_mask(tot);
      bitnames[std::string("AbsDrift")]         = bit_to_mask(absdrift);
      bitnames[std::string("DriftDt")]          = bit_to_mask(driftdt);
      bitnames[std::string("NullDriftVar")]     = bit_to_mask(nhdrift);
      bitnames[std::string("LongVal")]          = bit_to_mask(longval);
      bitnames[std::string("ANNProb")]          = bit_to_mask(annprob);
      bitnames[std::string("Added")]          = bit_to_mask(added);
      bitnames[std::string("GoodUDResid")]          = bit_to_mask(goodudresid);
      bitnames[std::string("GoodUTResid")]          = bit_to_mask(goodutresid);
      bitnames[std::string("GoodULResid")]          = bit_to_mask(goodulresid);
    }
    return bitnames;
  }
}
