#ifndef RecoDataProducts_StrawDigiFlag_hh
#define RecoDataProducts_StrawDigiFlag_hh
//
// Class to describe straw digi bits
//
// Original author David Brown
//
// Mu2e includes
#include "Offline/GeneralUtilities/inc/BitMap.hh"
#include <string>
#include <map>
#include <cstdint>
namespace mu2e {

  struct StrawDigiFlagDetail {
    typedef uint8_t mask_type;
    enum bit_type {energysel=0, dedxsel=1, corrupted=2, processed=7};
    // functions needed for the BitMap template
    static std::string const& typeName();
    static std::map<std::string,mask_type> const& bitNames();
    static mask_type bit_to_mask( bit_type b){ return 1<<b; }
 };
  typedef BitMap<StrawDigiFlagDetail> StrawDigiFlag;
  typedef std::vector<mu2e::StrawDigiFlag> StrawDigiFlagCollection;
}
#endif

