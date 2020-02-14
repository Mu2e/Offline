//
// Class to describe flag bits used for straw hits - D. Brown (LBL)
// 
#ifndef RecoDataProducts_BkgClusterFlag_hh
#define RecoDataProducts_BkgClusterFlag_hh

#include "GeneralUtilities/inc/BitMap.hh"
#include <string>
#include <map>

namespace mu2e 
{
  struct BkgClusterFlagDetail 
  {  
    // I need 32 bits for this class
    typedef unsigned mask_type;
    
    enum bit_type {stereo=0, update=1, unchanged=2, tight=4, loose=5, refined=10, bkg=15, iso=16};
    
    static std::string const& typeName();
    static std::map<std::string,mask_type> const& bitNames();
    static mask_type bit_to_mask(bit_type b) {return 1<<b;}
  };
  
  typedef BitMap<BkgClusterFlagDetail> BkgClusterFlag;

}
#endif

