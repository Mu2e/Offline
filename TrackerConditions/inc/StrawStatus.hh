#ifndef RecoDataProducts_StrawStatus_hh
#define RecoDataProducts_StrawStatus_hh
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
      absent=0, // physically absent
      nowire=1, // wire is removed or broken
      noHV=2, // wire is not attached to high voltage (blown fuse, ...)
      noLV=3, // no low voltage
      nogas=4, // gas flow is turned off or blocked or reduced
      noPreamp=5, // no signal into ADC or TDC
      noADC=6, // ADC not functioning
      noTDC=7, // TDC not functioning
      sparking=8, // straw generating microsparks
      noise=9, // straw electronics generating noise
      pickup=10, // straw noisy due to pickup
      suppress=11, // unspecified suppression
    };

// functions needed for the BitMap template
    static std::string const& typeName();
    static std::map<std::string,mask_type> const& bitNames();
// maximum track Id I can flag
    static unsigned _maxTrkId;
    static mask_type bit_to_mask( bit_type b){ return 1<<b; }
  };
  typedef BitMap<StrawStatusDetail> StrawStatus;
  typedef std::vector<mu2e::StrawStatus> StrawStatusCollection;
}
#endif

