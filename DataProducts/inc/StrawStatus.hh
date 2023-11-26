#ifndef DataProducts_StrawStatus_hh
#define DataProducts_StrawStatus_hh
//
// Class to describe flag bits used for defining straw (or panel or plane) status
//
// Original author David Brown (7/2020)
//
// Mu2e includes
#include "Offline/GeneralUtilities/inc/BitMap.hh"
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/StrawIdMask.hh"
#include <string>
#include <map>
namespace mu2e {

  struct StrawStatusDetail {
    // I need 32 bits for this class
    typedef unsigned mask_type;
    enum bit_type {
      absent=0, // straw (mylar + wire + gas) or panel or plane is physically absent
      nowire=1, // wire is removed or broken
      noHV=2, // Fuse blown (should be on both channels on this preamp) or HV problem
      noLV=3, // no low voltage (ie no electronics gain or signal)
      nogas=4, // gas flow is turned off
      lowgasgain=5, // gas gain low
      noHVPreamp=6, // HV preamp not installed or not connected to wire
      noCalPreamp=7, // Cal preamp not installed or not connected to wire
      noADC=8, // ADC not functioning
      noTDC=9, // TDC not functioning
      disabled=10, // readout of channel disabled in firmware
      sparking=11, // straw generating microsparks
      noise=12, // straw electronics generating noise
      pickup=13, // straw noisy due to pickup (ie adjacent noisy straw)
      suppress=14, // suppress signals from this straw for unspecified reasons
      nosurvey=15, // this element hasn't been surved which can affect its geometric accuracy
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

