#ifndef TrackerConditions_StrawStatus_hh
#define TrackerConditions_StrawStatus_hh
//
// Class to describe flag bits used for defining straw (or panel or plane) status
// 
// Original author David Brown (7/2020)
//
// Mu2e includes
#include "GeneralUtilities/inc/BitMap.hh"
#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/StrawIdMask.hh"
#include <string>
#include <map>
namespace mu2e {

struct StrawStatusDetail {
// I need 32 bits for this class
    typedef unsigned mask_type;
    enum bit_type {
      absent=0, // straw (mylar + wire + gas) or panel or plane is physically absent
      nowire=1, // wire is removed or broken
      noHV=2, // wire(s) not attached to high voltage (blown fuse, ...)
      noLV=3, // no low voltage (ie no electronics gain or signal)
      nogas=4, // gas flow is turned off
      lowgasgain=5, // gas gain low
      noPreamp=6, // no signal into ADC or TDC
      noADC=7, // ADC not functioning
      noTDC=8, // TDC not functioning
      sparking=9, // straw generating microsparks
      noise=10, // straw electronics generating noise
      pickup=11, // straw noisy due to pickup (ie adjacent noisy straw)
      suppress=12, // suppress signals from this straw for unspecified reasons
    };

// functions needed for the BitMap template
    static std::string const& typeName();
    static std::map<std::string,mask_type> const& bitNames();
    static mask_type bit_to_mask( bit_type b){ return 1<<b; }
  };
  typedef BitMap<StrawStatusDetail> StrawStatus;
  typedef std::vector<mu2e::StrawStatus> StrawStatusCollection;
}
#endif

