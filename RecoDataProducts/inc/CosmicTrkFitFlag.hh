#ifndef RecoDataProducts_CosmicTrkFitFlag_hh
#define RecoDataProducts_CosmicTrkFitFlag_hh
//Author: S Middleton
//Purpose: Class to describe flag bits for Cosmic Track Fitting

#include "GeneralUtilities/inc/BitMap.hh"
#include <string>
#include <map>
namespace mu2e {

  struct CosmicTrkFitFlagDetail {
    // I need 32 bits for this class
    typedef unsigned mask_type;
    // The first 16 describe various success conditions, the last 16 various failure modes
    enum bit_type {hitsOK=0,seedOK,seedConverged, StraightTrackOK,StraightTrackConverged, StraightTrackInit};

    // functions needed for the BitMap template
    static std::string const& typeName();
    static std::map<std::string,mask_type> const& bitNames();
    static mask_type bit_to_mask( bit_type b){ return 1<<b; }
  };
  typedef BitMap<CosmicTrkFitFlagDetail> CosmicTrkFitFlag;

}
#endif

