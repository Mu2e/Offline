#ifndef RecoDataProducts_StrawFlag_hh
#define RecoDataProducts_StrawFlag_hh
//
// Class to describe flag bits used for straw materials
//
// Original author David Brown 1/8/2025
//
// Mu2e includes
#include "Offline/GeneralUtilities/inc/BitMap.hh"
#include <string>
#include <map>
namespace mu2e {

  struct StrawFlagDetail {
// I need 32 bits for this class
    typedef unsigned mask_type;
    enum bit_type {active=0, hashit=1, activehit=2, drifthit=3};
    // functions needed for the BitMap template
    static std::string const& typeName();
    static std::map<std::string,mask_type> const& bitNames();
// maximum track Id I can flag
    static unsigned _maxTrkId;
    static mask_type bit_to_mask( bit_type b){ return 1<<b; }
  };
  typedef BitMap<StrawFlagDetail> StrawFlag;
  typedef std::vector<mu2e::StrawFlag> StrawFlagCollection;
}
#endif
