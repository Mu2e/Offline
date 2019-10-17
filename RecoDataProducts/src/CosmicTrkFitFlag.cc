//Author: S Middleton
//Purpose: Flags for Cosmic Track fitting
#include <stdexcept>
#include <iostream>
#include <stdio.h>
#include "RecoDataProducts/inc/CosmicTrkFitFlag.hh"
namespace mu2e {

  std::string const& CosmicTrkFitFlagDetail::typeName() {
    static std::string type("CosmicTrkFitFlag");
    return type;
  }

  std::map<std::string,CosmicTrkFitFlagDetail::mask_type> const& CosmicTrkFitFlagDetail::bitNames() {
    static std::map<std::string,mask_type> bitnames;
    if(bitnames.size()==0){
      
      bitnames[std::string("HitsOK")]             = bit_to_mask(hitsOK);
      bitnames[std::string("SeedOK")]              = bit_to_mask(seedOK);
      bitnames[std::string("StraightTrackOK")]        = bit_to_mask(StraightTrackOK);
      bitnames[std::string("StraightTrackConverged")]   = bit_to_mask(StraightTrackConverged);
      bitnames[std::string("StraightTrackInit")]     = bit_to_mask(StraightTrackInit);

    }
    return bitnames;
  }

}
