#ifndef RecoDataProducts_WHSMask_hh
#define RecoDataProducts_WHSMask_hh
//
// Class to describe a mask of WireHitState states
// Original author David Brown
//
// Mu2e includes
#include "Offline/GeneralUtilities/inc/BitMap.hh"
#include <string>
#include <map>
namespace mu2e {

  struct WHSMaskDetail {
    typedef unsigned mask_type;
    enum bit_type {inactive=0,null,drift};
    // functions needed for the BitMap template
    static std::string const& typeName();
    static std::map<std::string,mask_type> const& bitNames();
    static mask_type bit_to_mask( bit_type b){ return 1<<b; }
  };
  typedef BitMap<WHSMaskDetail> WHSMask;
}
#endif
