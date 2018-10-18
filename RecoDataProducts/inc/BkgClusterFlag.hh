#ifndef RecoDataProducts_BkgClusterFlag_hh
#define RecoDataProducts_BkgClusterFlag_hh
//
// Class to describe flag bits used for straw hits
// 
// $Id: BkgClusterFlag.hh,v 1.4 2013/04/04 01:08:20 brownd Exp $
// $Author: brownd $
// $Date: 2013/04/04 01:08:20 $
//
// Original author David Brown
//
// Mu2e includes
#include "GeneralUtilities/inc/BitMap.hh"
#include <string>
#include <map>
namespace mu2e {

  struct BkgClusterFlagDetail {
    // I need 32 bits for this class
    typedef unsigned mask_type;
    // different classification bits
    enum bit_type {stereo=0, tight=4, loose=5, refined=10, bkg=15, iso=16};
    // functions needed for the BitMap template
    static std::string const& typeName();
    static std::map<std::string,mask_type> const& bitNames();
    static mask_type bit_to_mask( bit_type b){ return 1<<b; }
  };
  typedef BitMap<BkgClusterFlagDetail> BkgClusterFlag;

}
#endif

